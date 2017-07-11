#!/usr/bin/env python
import sys
import os
import pysam
from lib.read_stuff import get_coordinate
from lib.file_stuff import get_bamfile, get_output

def no_recal_tag(read):
    return not (read.has_tag('BD') or read.has_tag('BI'))

def usage():

    msg = "Usage: {} input.bam output.bam [dups.bam]".format(sys.argv[0])
    msg += '''

Removes unmapped reads without GATK base quality 
score recalibration tags (BD and BI) where the 
mate is mapped (pairs where both reads are 
unmapped are kept). This is an attempt to 
mitigate a problem with duplicate unmapped reads
being written to output files by some bcbio 
pipelines. 

Expects input alignment file as the first 
argument, output bam filename as the second. 
Optionally, a filename for filtered reads can be 
given as the third argument.

The functionality of this script is superceded by
the identify_dups.py script, but this script may
be useful for checking/testing.
'''
    sys.exit(msg)

if __name__ == '__main__':
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        usage()
    bam = sys.argv[1]
    out = sys.argv[2]
    dups = None
    if len(sys.argv) > 3:
        dups = sys.argv[3]
    for fn in (out, dups):
        if fn is not None and os.path.exists(fn):
            sys.exit("Output file '{}' exists - please delete or choose"
                     .format(fn) + " another name.")
    bamfile = get_bamfile(bam)
    outfile = get_output(out, bamfile)
    if dups is not None:
        dupfile = get_output(dups, bamfile)
    n = 0
    skipped = 0
    prog_string = ''
    for read in bamfile.fetch(until_eof=True):
        n += 1
        if (read.is_unmapped and no_recal_tag(read) and 
            not read.mate_is_unmapped): 
            skipped += 1
            if dupfile is not None:
                dupfile.write(read)
        else:
           outfile.write(read)
        if not n % 10000:
            coord = get_coordinate(read)
            sys.stderr.write("\r" + " " * len(prog_string))
            prog_string = "{:,} read, {:,} filtered,  at {}".format(n, skipped, 
                                                                    coord)
            sys.stderr.write("\r" + prog_string)
    sys.stderr.write(
                  "\rFinished filtering: {:,} alignments read, {:,} filtered\n"
                  .format(n, skipped))
    bamfile.close()
    outfile.close()
    dupfile.close()
