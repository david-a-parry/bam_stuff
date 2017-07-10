#!/usr/bin/env python
import sys
import os
import pysam
from collections import defaultdict,OrderedDict
import argparse

def get_parser():
    '''Get ArgumentParser'''
    parser = argparse.ArgumentParser(
                  description='''   Identify and filter erroneously duplicated 
                                    reads from a BAM. Keeps only the first 
                                    encountered read where more than one share 
                                    the same read ID and bitwise FLAG.''')
    required_args = parser.add_argument_group('Required Arguments')
    optional_args = parser.add_argument_group('Optional Arguments')
    required_args.add_argument('-i', '--bam', '--input', required=True, 
                               metavar='IN.bam', help='''Input BAM filename''')
    required_args.add_argument('-o', '--output', metavar='OUT.bam', 
                               required=True,
                               help='''Write cleaned output BAM to this file. 
                                    ''')
    optional_args.add_argument('-d', '--dup_bam', metavar='DUPS.bam', 
                               help='''Optional BAM file for duplicated reads.
                               ''')
    optional_args.add_argument('-r', '--ref', metavar='REF', 
                               help='''Reference fasta file. Required if using
                                       CRAM input.''')
    optional_args.add_argument('-b', '--buffer_size', metavar='NUM_READS', 
                               type=int, default=500, 
                               help='''Number of reads to store in memory to 
                                       check as duplicates against the current 
                                       read. Only the last N reads will be held
                                       for comparison. Too few may mean that
                                       duplicates are not detected while a 
                                       large number will result in slower 
                                       runtime and greater memory usage.
                                       Default=500.''')
    return parser
    
 
def get_coordinate(read):
    if read.reference_id >= 0:
        coord = "{}:{:,}".format(read.reference_name, read.reference_start)
    else:
        coord = "*/*"
    return coord

def get_duplicates(read, cache):
    dups = []
    if read.query_name in cache:
        for cr in cache[read.query_name]:
            if is_dup(read, cr):
                dups.append(cr)
    return dups

def is_dup(read, other):
    # cache key is query_name so we know read.query_name == other.query_name
    return (read.flag == other.flag)
    #return (read.flag == other.flag and 
    #    read.reference_start == other.reference_start and 
    #    read.seq == other.seq)
        # not quite exhaustive, but surely have to be duplicates if 
        # all these are true. Maybe just name and Flag is enough?

def filter_dups(bam, output, dup_bam=None, ref=None, buffer_size=500):
    for fn in (output, dup_bam):
        if os.path.exists(fn):
            sys.exit("Output file '{}' exists - please delete or choose"
                     .format(fn) + " another name.")
    if buffer_size < 1:
        sys.exit("--buffer_size must be greater than 0")
    bmode = 'rb'
    if bam.endswith(('.sam', '.SAM')):
        bmode = 'r'
    elif bam.endswith(('.cram', '.CRAM')):
        bmode = 'rc'
    bamfile = pysam.AlignmentFile(bam, bmode)
    outfile = pysam.AlignmentFile(output, "wb", template=bamfile)
    dupfile = None
    if dup_bam is not None:
        dupfile = pysam.AlignmentFile(dup_bam, "wb", template=bamfile)
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
        cache_keys[read.query_name] = 1
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

