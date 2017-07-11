#!/usr/bin/env python
import sys
import os
import pysam
from collections import defaultdict,OrderedDict
import argparse
from lib.read_stuff import get_coordinate
from lib.file_stuff import get_bamfile, get_output

def get_parser():
    '''Get ArgumentParser'''
    parser = argparse.ArgumentParser(
                  description='''   Identify and filter erroneously duplicated 
                                    reads from a BAM. Keeps only the first 
                                    encountered read where more than one share 
                                    the same read ID and bitwise FLAG.''')
    required_args = parser.add_argument_group('Required Arguments')
    optional_args = parser.add_argument_group('Optional Arguments')
    required_args.add_argument('bam', metavar='IN', 
                               help=''' Input BAM filename''')
    optional_args.add_argument('-o', '--output', metavar='OUT', 
                               help=''' Write cleaned output BAM to this file. 
                                        Default is to write SAM format to 
                                        STDOUT''')
    optional_args.add_argument('-d', '--dup_bam', metavar='DUPS.bam', 
                               help=''' Optional BAM file for duplicated reads.
                               ''')
    optional_args.add_argument('-r', '--ref', metavar='REF', 
                               help=''' Reference fasta file. Required if using
                                        CRAM input.''')
    optional_args.add_argument('-b', '--buffer_size', metavar='NUM_READS', 
                               type=int, default=500, 
                               help=''' Number of reads to store in memory to 
                                        check as duplicates against the current 
                                        read. Only the last N reads will be 
                                        held for comparison. Too few may mean 
                                        that duplicates are not detected while 
                                        large number may result in slower 
                                        runtime and greater memory usage. The 
                                        ideal number will depend on depth of 
                                        coverage (i.e. how many reads may 
                                        have the same start coordinate), with 
                                        the default setting being adequate for
                                        30X coverage. Default=500.''')
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

def filter_dups(bam, output=None, dup_bam=None, ref=None, buffer_size=500):
    for fn in (output, dup_bam):
        if fn is not None and os.path.exists(fn):
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
    for read in bamfile.fetch(until_eof=True):
        n += 1
        dups = get_duplicates(read, cache)
        if dups:
            d += 1
            if dupfile is not None:
                for dr in dups:
                    dupfile.write(dr)
                dupfile.write(read)
        else:
            outfile.write(read)
            cache[read.query_name].append(read)
            cache_keys[read.query_name] = None
            if len(cache_keys) > buffer_size:
                first = list(cache_keys)[0]
                del cache[first]
                del cache_keys[first]
        if not n % 10000:
            coord = get_coordinate(read)
            sys.stderr.write("\r" + " " * len(prog_string))
            prog_string = ("{:,} read, {:,} duplicates identified".format(n, d)  
                          + " at {}".format(coord))
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

