#!/usr/bin/env python
import sys
import os
import pysam
from lib.read_stuff import get_coordinate
from lib.file_stuff import get_bamfile, get_output

def no_recal_tag(read):
    return not (read.has_tag('BD') or read.has_tag('BI'))

def usage():
        msg = "Usage: {} input.bam [ref.fa]" .format(sys.argv[0])
        msg += '''

Outputs number of mapped and unmapped reads with 
and without the GATK base quality recalibration
tags 'BD' or 'BI'.

First argument must be the input BAM/SAM/CRAM 
file. The optional second argument is a reference
fasta file, which is required if using CRAM input.

The purpose of this script is to investigate BAM
files that have gone through bcbio pipelines which
result in not all reads having the expected BQSR
tags.

'''
        sys.exit(msg)

if __name__ == '__main__':
    if len(sys.argv) < 2 or len(sys.argv) > 3: 
        usage()
    bam = sys.argv[1]
    ref = None
    if len(sys.argv) > 2:
        ref = sys.argv[2]
    filtered_reads = None
    bamfile = get_bamfile(bam)
    n = 0
    supp = 0
    prog_string = ''
    categories = {'with_tag': {'mapped': 0, 'unmapped': 0}, 
                  'no_tag': {'mapped': 0, 'unmapped': 0}}
    for read in bamfile.fetch(until_eof=True, reference=ref):
        n += 1
        k1 = 'with_tag'
        k2 = 'mapped'
        if no_recal_tag(read):
            k1 = 'no_tag'
        if read.is_unmapped:
            k2 = 'unmapped'
        categories[k1][k2] += 1
        if not n % 10000:
            coord = get_coordinate(read)
            sys.stderr.write("\r" + " " * len(prog_string))
            prog_string = ("{:,} read, ".format(n) + "at {}".format(coord))
            sys.stderr.write("\r" + prog_string)
    sys.stderr.write("\rFinished: {:,} alignments read\n".format(n))
    bamfile.close()
    print(str.join("\t", ("Type", "With_Recal_Tag", "No_Recal_Tag")))
    print(str.join("\t", ('mapped',
                   str(categories['with_tag']['mapped']),
                   str(categories['no_tag']['mapped']))))
    print(str.join("\t", ('unmapped',
                   str(categories['with_tag']['unmapped']),
                   str(categories['no_tag']['unmapped']))))
