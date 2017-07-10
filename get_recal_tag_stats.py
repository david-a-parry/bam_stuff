#!/usr/bin/env python
import sys
import os
import pysam

def no_recal_tag(read):
    return not (read.has_tag('BD') or read.has_tag('BI'))

def get_coordinate(read):
    if read.reference_id >= 0:
        coord = "{}:{:,}".format(read.reference_name, read.reference_start)
    else:
        coord = "*/*"
    return coord

if __name__ == '__main__':
    if len(sys.argv) < 2 or len(sys.argv) > 3: 
        sys.exit("Usage: {} input.bam [ref.fa]" .format(sys.argv[0]))
    bam = sys.argv[1]
    ref = None
    if len(sys.argv) > 2:
        ref = sys.argv[2]
    filtered_reads = None
    bmode = 'rb'
    if bam.endswith(('.sam', '.SAM')):
        bmode = 'r'
    elif bam.endswith(('.cram', '.CRAM')):
        bmode = 'rc'
    bamfile = pysam.AlignmentFile(bam, bmode)
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
