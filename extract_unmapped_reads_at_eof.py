#!/usr/bin/env python
import sys
import os
import pysam
import argparse
from lib.read_stuff import get_coordinate
from lib.file_stuff import get_bamfile, get_output, seek_back_til_reads

def get_parser():
    '''Get ArgumentParser'''
    parser = argparse.ArgumentParser(
       formatter_class=argparse.RawTextHelpFormatter,
        description=''' 
Extract unmapped reads from the end of a BAM file.''',
        epilog='''  

    Read all unmapped reads (without a mapped mate) from a BAM file. 

    While samtools will extract all reads for a given region/contig it
    does not provide a way to seek to unmapped reads at the end of a 
    file. This script provides a way to retrieve all unmapped reads from 
    the end of a BAM file. Note, that if your reads are paired, this 
    will not retrieve unmapped reads where the mate is mapped. That can
    instead be acheived simply using samtools (for example:
                                               samtools -f 4 input.bam) 

''' )
    required_args = parser.add_argument_group('Required Arguments')
    optional_args = parser.add_argument_group('Optional Arguments')
    required_args.add_argument('-i', '--bam', '--input', required=True, 
                               metavar='IN', 
                               help='''Input BAM filename''')
    optional_args.add_argument('-o', '--output', metavar='OUT', 
                               help=
'''Write cleaned output BAM to this file. Default
behaviour is to write SAM format to STDOUT.

''')
    optional_args.add_argument('-r', '--ref', metavar='REF', 
                               help=
'''Reference fasta file. Required if using CRAM input.''')
    return parser

def extract(bam, output=None, ref=None):
    if output is not None and os.path.exists(output):
        sys.exit("Output file '{}' exists - please delete or choose"
                 .format(output) + " another name.")
    bamfile = get_bamfile(bam)
    outfile = get_output(output, bamfile)
    n = 0
    u = 0
    prog_string = ''
    seek_back_til_reads(bamfile, ref=ref)
    sys.stderr.write("File position should now be at unmapped pairs\n")
    for read in bamfile.fetch(until_eof=True):
        u += 1
        outfile.write(read)
        if not u % 1000:
            sys.stderr.write("\r{:,} unmapped reads written".format(u))
    sys.stderr.write("\n{:,} unmapped reads with unmapped mates found\n"
                     .format(u))
    bamfile.close()
    outfile.close()
 
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    extract(**vars(args))

