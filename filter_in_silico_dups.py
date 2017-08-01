#!/usr/bin/env python
import sys
import os
import pysam
from collections import defaultdict,OrderedDict
import argparse
from lib.read_stuff import get_coordinate, get_rname
from lib.file_stuff import get_bamfile, get_output

def get_parser():
    '''Get ArgumentParser'''
    parser = argparse.ArgumentParser(
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    usage='%(prog)s INPUT [options]',
                    description=''' 
Identify and filter erroneously duplicated reads from a BAM. Keeps only the 
first encountered read where more than one share the same read ID and bitwise 
FLAG.''',

                    epilog='''

A note on the program logic:
    
In theory, for all primary alignments the read IDs and bitwise flags should be
unique, meaning we can filter duplicate reads using these features. For
secondary or supplementary alignments bitwise flags need not be unique, so these
alignments are only considered duplicates if POS and CIGAR fields are identical.

Note, that for primary alignments, we cannnot consider unmapped reads with
different POS values as necessarily distinct as this value may change if the
mapped mate has undergone indel realignment. This means an arbitrary sized read
buffer is necessary rather than just retaining reads mapped to the same
coordinate in memory.''')

    required_args = parser.add_argument_group('Required Arguments')
    optional_args = parser.add_argument_group('Optional Arguments')
    required_args.add_argument('bam', metavar='INPUT', 
                               help=''' Input BAM filename''')
    optional_args.add_argument('-o', '--output', metavar='OUT', 
                               help=''' Write cleaned output BAM to this file. 
                                        Default is to write SAM format to 
                                        STDOUT''')
    optional_args.add_argument('-d', '--dup_bam', metavar='DUPS', 
                               help=''' Optional BAM file for duplicated reads.
                               ''')
    optional_args.add_argument('-r', '--ref', metavar='REF', 
                               help=''' Reference fasta file. May be required   
                                        if using CRAM input.''')
    optional_args.add_argument('-w', '--window_size', metavar='DISTANCE',
                               type=int, default=100, 
                               help=''' Check for duplicate reads with aligned
                                        coordinates this far apart or closer. 
                                        For reads where both pairs are unmapped
                                        --buffer_size is used to determine the
                                        number of reads to keep instead. This 
                                        is meant to allow for potential 
                                        adjustments in the aligned coordinate 
                                        due to GATK indel realignment while 
                                        preserving memory usage.
                                        Default=100''')
    optional_args.add_argument('-b', '--buffer_size', metavar='NUM_READS', 
                               type=int, default=100, 
                               help=''' Number of reads to store in memory to 
                                        check as duplicates against the current 
                                        read. Only the last N reads will be 
                                        held for comparison if the first and 
                                        last span a greater distance than 
                                        --window_size or else if processing
                                        unmapped pairs. Default=100.''')
    optional_args.add_argument('-f', '--force', action='store_true', 
                               help=''' Overwite existing output files.''')
    return parser
    
def get_duplicates(read, cache):
    dups = []
    if read.query_name in cache:
        for cr in cache[read.query_name]:
            if is_dup(read, cr):
                dups.append(cr)
    return dups

def is_dup(read, other):
    # cache key is query_name so we know read.query_name == other.query_name
    if read.flag == other.flag:
        if read.is_secondary or read.is_supplementary:
            return (read.reference_start == other.reference_start and 
                    read.cigarstring == other.cigarstring)
        else:
            return True
    return False

def check_sorted(read, prev, seen_contigs):
    if prev is None:
        return
    if read.reference_id == prev.reference_id:
        if read.reference_start < prev.reference_start:
            sys.exit("ERROR: Unsorted positions: {} after {}"
                     .format(get_coordinate(read), get_coordinate(prev)))
    else:
        if read.reference_id in seen_contigs:
            sys.exit("ERROR: Unsorted positions: encountered contig " + 
                     " {} before and after {}" 
                     .format(get_rname(read), get_rname(prev)))
        seen_contigs.add(get_rname(read))

def pop_cache(read, prev, window, size, cache, cache_keys):
    ''' 
        If we are at a new position and read name is not in cache 
        already, remove the first N cache entries to reduce cache size
        to buffer_size.
    '''
    to_remove = []
    if (read.reference_start != prev.reference_start or
        read.reference_id != prev.reference_id or read.reference_id == -1):
        if read.query_name not in cache_keys:
            r = len(cache_keys) - size
            if r > 0:
                ck = list(cache_keys)
                check_window = False
                if (read.reference_id == prev.reference_id and 
                    read.reference_id != -1):
                    #keep mapped reads close to each other
                    check_window = True
                    window_end = cache[ck[-1]][-1].reference_start
                    window_ref_id = cache[ck[-1]][-1].reference_id
                for k in ck[0:r]:
                    if check_window:
                        if (cache[k][-1].reference_id == window_ref_id and
                           window_end - cache[k][-1].reference_start < window):
                            #bail out if first read is within window bp of last
                            break
                    del cache[k]
                    del cache_keys[k]

def filter_dups(bam, output=None, dup_bam=None, ref=None, window_size=100,
                buffer_size=100, force=False):
    for fn in (output, dup_bam):
        if fn is not None and os.path.exists(fn) and not force:
            sys.exit("Output file '{}' exists - please delete or choose"
                     .format(fn) + " another name.")
    if buffer_size < 1:
        sys.exit("--buffer_size must be greater than 0")
    bamfile = get_bamfile(bam)
    outfile = get_output(output, bamfile) 
    dupfile = None
    if dup_bam is not None:
        dupfile = get_output(dup_bam, bamfile)
    n = 0
    d = 0
    prog_string = ''
    cache = defaultdict(list)
    cache_keys = OrderedDict()
    contigs = set()
    buffer_full = False
    prev_read = None
    for read in bamfile.fetch(until_eof=True):
        n += 1
        check_sorted(read, prev_read, contigs)
        dups = get_duplicates(read, cache)
        if dups:
            d += 1
            if dupfile is not None:
                for dr in dups:
                    dupfile.write(dr)
                dupfile.write(read)
        else:
            outfile.write(read)
            if buffer_full:
                pop_cache(read, prev_read, window_size, buffer_size, cache, 
                          cache_keys)
            elif len(cache_keys) >= buffer_size:
                buffer_full = True
                pop_cache(read, prev_read, window_size, buffer_size, cache, 
                          cache_keys)
            cache[read.query_name].append(read)
            cache_keys[read.query_name] = None
        prev_read = read
        if not n % 10000:
            sys.stderr.write("\r" + " " * len(prog_string))
            prog_string = ("{:,} read, {:,} duplicates identified".format(n, d)  
                          + " at {}".format(get_coordinate(read)))
            sys.stderr.write("\r" + prog_string)
    sys.stderr.write("\nFinished: {:,} alignments read, {:,} filtered\n"
                     .format(n, d))
    bamfile.close()
    outfile.close()
    if dupfile is not None:
        dupfile.close()

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    filter_dups(**vars(args))

