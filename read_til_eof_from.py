#!/usr/bin/env python
import sys
import os
import pysam

def get_rname(read):
    if read.reference_id >= 0:
        return read.reference_name
    return '*'

def usage():
    msg = "\n    Usage: {} input.bam contig [ref.fa]" .format(sys.argv[0])
    msg += '''

    Read from a given contig to the end of an indexed
    BAM/CRAM file, including unmapped reads. Output is
    written to STDOUT in SAM format. A reference is 
    required as the third argument if CRAM format is 
    used.

    While samtools will extract all reads for a given
    region/contig it does not provide a way to seek to
    unmapped reads at the end of a file. This script 
    provides a way to retrieve ALL reads from a given 
    contig onwards, which may be useful in certain 
    circumstances (mostly debugging or exploratory 
    purposes).

    If you are ONLY interested in unmapped reads, see
    the extract_unmapped_pairs.py script instead.

'''
    sys.exit(msg)

if __name__ == '__main__':
    if (len(sys.argv) < 3 or len(sys.argv) > 4):
        usage()
    bam = sys.argv[1]
    contig = sys.argv[2]
    ref = None
    if len(sys.argv) > 3:
        ref = sys.argv[3]
    bmode = 'rb'
    if bam.endswith(('.sam', '.SAM')):
        bmode = 'r'
    elif bam.endswith(('.cram', '.CRAM')):
        bmode = 'rc'
    bamfile = pysam.AlignmentFile(bam, bmode)
    tid = bamfile.get_tid(contig)
    if tid < 0:
        sys.exit("Could not find contig '{}' in bam".format(contig))
    outfile = pysam.AlignmentFile('-', "w", template=bamfile)
    sys.stderr.write("Seeking to reference ({})"
                     .format(bamfile.get_reference_name(tid)))
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
 
