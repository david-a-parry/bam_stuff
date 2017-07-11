#!/usr/bin/env python
import sys
import os
import pysam
from lib.read_stuff import get_rname
from lib.file_stuff import get_bamfile, get_output
import argparse

def get_parser():
    '''Get ArgumentParser'''
    parser = argparse.ArgumentParser(
            description=''' Read from a given contig to the end of an indexed
                            BAM/CRAM file, including unmapped reads.''',
            epilog='''  While samtools will extract all reads for a given
                        region/contig it does not provide a way to seek to
                        unmapped reads at the end of a file. This script 
                        provides a way to retrieve ALL reads from a given 
                        contig onwards, which may be useful in certain 
                        circumstances (mostly debugging or exploratory 
                        purposes).

                        If you are ONLY interested in unmapped reads, see
                        the extract_unmapped_pairs.py script instead.''')

    required_args = parser.add_argument_group('Required Arguments')
    optional_args = parser.add_argument_group('Optional Arguments')
    required_args.add_argument('bam', metavar='BAM', 
                               help=''' Input BAM filename''')
    required_args.add_argument('contig', metavar='CONTIG', 
                               help=''' Name of contig to start retrieving 
                                        reads from''')
    optional_args.add_argument('-o', '--output', metavar='OUT', 
                               help=''' Write output to this file. Default 
                                        behaviour is to write SAM format to 
                                        STDOUT.''')
    optional_args.add_argument('-r', '--ref', metavar='REF', 
                               help=''' Reference fasta file. Required if using 
                                        CRAM input.''')
    return parser


def read_til_eof(bam, contig, output=None, ref=None):
    bamfile = get_bamfile(bam)
    tid = bamfile.get_tid(contig)
    if tid < 0:
        sys.exit("Could not find contig '{}' in bam".format(contig))
    outfile = get_output(output, bamfile)
    sys.stderr.write("Seeking to reference {}"
                     .format(bamfile.get_reference_name(tid)))
    sys.stderr.write("\nAt contig {}".format(bamfile.get_reference_name(tid)))
    for read in bamfile.fetch(tid=tid, reference=ref):
        outfile.write(read)
    #get remaining reads
    prev_rname = contig
    for read in bamfile.fetch(until_eof=True, reference=ref):
        rname = get_rname(read)
        if rname != prev_rname:
            sys.stderr.write("\nAt contig " + rname)
            prev_rname = rname
        outfile.write(read)
    sys.stderr.write("\nFinished\n")
    bamfile.close()
    outfile.close()
 
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    read_til_eof(**vars(args))
