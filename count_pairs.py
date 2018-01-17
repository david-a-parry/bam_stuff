#!/usr/bin/env python3
import sys
import os
import pysam
import argparse
import logging
from lib.read_stuff import get_coordinate, get_rname
from lib.file_stuff import get_bamfile, get_output

logger = logging.getLogger("CountPairs")
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
                '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler()
ch.setLevel(logger.level)
ch.setFormatter(formatter)
logger.addHandler(ch)

def get_parser():
    '''Get ArgumentParser'''
    parser = argparse.ArgumentParser(
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    usage='%(prog)s BAM [options]',
                    description=''' 
Count the number of read pairs where at least one read of the pair satisfies 
the given criteria.''',)

    parser.add_argument('bam', metavar='BAM', 
                               help=''' Input BAM file.''')
    parser.add_argument('-m', '--mapping_quality', metavar='MAPQ', type=int, 
                        default=30,
                               help=''' Minimum mapping quality of reads to
                                        consider. Default=30.''')
    parser.add_argument('-c', '--collated', 
                               help=''' Prefix for temporary collated BAM file.
                                        Defaults to input BAM name plus 
                                        ".collated"''')
    parser.add_argument('-t', '--threads', type=int, default=1,
                               help=''' Number of threads to use for collating
                                        reads. Default=1''')
    parser.add_argument('-d', '--duplicates', action='store_true',
                               help=''' Count alignments marked as duplicates.
                                        Default behaviour is to ignore 
                                        duplicates.''')
    parser.add_argument('-s', '--secondary', action='store_true',
                               help=''' Count secondary/supplementary 
                                        alignments as well as primary 
                                        alignments.''')
    return parser

def count_pairs(bam, mapping_quality=30, collated=None, threads=1, 
                duplicates=False, secondary=False):
    if collated is None:
        collated = bam + ".collated"
    if os.path.exists(collated + '.bam'):
       sys.exit("Error: collated output '{}.bam' already exists".format(
                collated))
    tmp_bam = collate_bam(bam, collated, threads)
    last_counted = ''
    p = 0
    for read in tmp_bam.fetch(until_eof=True):
        if read.query_name == last_counted: 
            continue
        if not secondary and (read.is_secondary or read.is_supplementary):
            continue
        if not duplicates and read.is_duplicate: 
            continue
        if read.mapping_quality >= mapping_quality:
            last_counted = read.query_name
            p += 1
    print(p)
    tmp_bam.close()
    os.remove(collated + '.bam')

def collate_bam(bam, output, threads):
    logger.info("Collating cached reads by read ID...")
    args = ['--output-fmt', 'BAM', bam, output]
    if threads > 1:
        args += ['-@', str(threads - 1)]
    pysam.collate(*args)
    logger.info("Finished collating cached reads.")
    return pysam.AlignmentFile(output + '.bam', 'rb')
    

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    count_pairs(**vars(args))

