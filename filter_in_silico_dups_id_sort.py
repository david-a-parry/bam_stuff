#!/usr/bin/env python
import sys
import os
from subprocess import PIPE, Popen
import argparse
from collections import defaultdict
import re
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
mapped mate has undergone indel realignment.''')

    parser.add_argument('bam', metavar='INPUT', 
                               help=''' Input BAM filename. Input MUST be 
                                        sorted by read name.''')
    parser.add_argument('-o', '--output', metavar='OUT', required=True,
                               help=''' Write cleaned output BAM to this file. 
                                        Default is to write SAM format to 
                                        STDOUT''')
    parser.add_argument('-d', '--dup_bam', metavar='DUPS', 
                               help=''' Optional BAM file for duplicated reads.
                               ''')
   # parser.add_argument('-r', '--ref', metavar='REF', 
   #                            help=''' Reference fasta file. May be required   
   #                                     if using CRAM input.''')
    parser.add_argument('-f', '--force', action='store_true', 
                               help=''' Overwite existing output files.''')
    return parser
    
def is_dup(read, other):
    if read.flag == other.flag:
        if read.is_secondary or read.is_supplementary:
            return (read.reference_start == other.reference_start and 
                    read.cigarstring == other.cigarstring)
        else:
            return True
    return False

def check_dups(cache):
    '''
        Identify duplicate reads in cache. 

        Returns two lists of reads. The first is a list of non-duplicated reads,
        the second is a list of filtered duplicates.
    '''
    non_dups = []
    dups = []
    read1s = []
    read2s = []
    singletons = []
    ndups = 0
    for r in cache:
        if r.is_paired:
            if r.is_read1:       #read1
                read1s.append(r)
            elif r.is_read2:    #read2
                read2s.append(r)
            else:
                sys.exit("ERROR: Malformed BAM - read is paired but " + 
                         "is neither set as first or second in pair for " +
                         "read {}".format(read.tostring(bamfile)))
        else:
            singletons.append(r)
    if singletons:
        singles, empty, supp = assign_primary_and_supplementary(singletons)
        ndups += add_non_dups_and_dups(singletons, non_dups, dups)
    r1s, r2s, psupp = assign_primary_and_supplementary(read1s, read2s)
    if psupp:
        ndups += add_non_dups_and_dups(psupp, non_dups, dups)
    if len(r1s) > 1:
        dups.extend(r1s)
        ndups += 1
    if len(r2s) > 1:
        dups.extend(r2s)
        ndups += 1
    if len(r1s) > 0 and len(r2s) > 0:
        non_dups.extend(get_non_dup_pair(r1s, r2s))
    elif len(r1s) > 0:
        non_dups.append(choose_best_read(r1s))
    elif len(r2s) > 0:
        non_dups.append(choose_best_read(r2s))
    return non_dups, dups, ndups


def get_non_dup_pair(r1s, r2s):
    if len(r1s) == 1 and len(r2s) == 1:
        return (r1s[0], r2s[0])
    pairs = []
    for r1 in r1s:
        for r2 in r2s:
            if (r1.flag & 4 and not r2.flag & 4 or  # unmapped and mapped 
                r2.flag & 4 and not r1.flag & 4):   # - expect same coordinate,
                #coords may differ if mapped read underwent indel realignment
                #and unmapped read originated from before indel realignment
                if (r1.reference_start == r2.reference_start and 
                    r1.reference_id == r2.reference_id):
                    pairs.append((r1, r2))
            elif not r1.flag & 4 and not r2.flag & 4: #both mapped 
                #check mate coords match
                if (r1.reference_start == r2.next_reference_start and 
                    r1.reference_id == r2.next_reference_id):
                    pairs.append((r1, r2))
            else:#both unmapped
                pairs.append((r1, r2))
    best_pair = None
    for p in pairs:
        if has_recal_tag(p[0]) and has_recal_tag(p[1]):
            best_pair = (p[0], p[1])
            break
        elif has_oc_tag(p[0]) or has_oc_tag(p[1]):
            best_pair = (p[0], p[1])
    if best_pair is None and pairs:
        best_pair = (pairs[0][0], paris[0][1])
    return best_pair

def add_non_dups_and_dups(reads, non_dups, dups):
    sdups,snon_dups = get_dup_reads(reads)
    non_dups.extend(snon_dups)
    if sdups:
        non_dups.extend(choose_best_dup(sdups))
        dups.extend(sdups)
    return len(sdups)

def assign_primary_and_supplementary(read1s, read2s=[]):
    pairs = defaultdict(dict)
        #dict of read 1 flag to dict of read1 and read2 lists
    suppl = [] #supplementary/non-primary alignments
    pr1 = [] #primary alignments for read1
    pr2 = []
    for r1 in read1s:
        if r1.flag & 256 or r1.flag & 2048: #not primary/supplementary 
            suppl.append(r1)
        else:
            pr1.append(r1)
    for r2 in read2s:
        if r2.flag & 256 or r2.flag & 2048: #not primary/supplementary 
            suppl.append(r2)
        else:
            pr2.append(r2)
    return pr1, pr2, suppl

def get_dup_reads(reads):
    ''' 
        Returns a list of non duplicate read (binary) as first item and
        a dict of flags to lists of ReadSummary objects as the second
    '''
    dup_sets = defaultdict(list)
    non_dups = []
    for i in range(len(reads)):
        is_dup = False
        for j in range(len(reads)):
            if i == j: continue
            if reads[i].flag == reads[j].flag:
                if reads[i].flag & 256 or reads[i].flag & 2048: 
                    #not primary/supplementary 
                    if (reads[i].reference_id == reads[j].reference_id and 
                        reads[i].reference_start == reads[j].reference_start):
                        #same coord
                        is_dup = True
                        dup_sets[reads[i].flag].append(reads[i])
                else:       
                    is_dup = True
                    dup_sets[reads[i].flag].append(reads[i])
            if is_dup:
                break
        if not is_dup:
            non_dups.append(reads[i])
    return dup_sets, non_dups

def choose_best_dup(dups):
    best_dups = []
    for k,v in dups.items(): #key is bitwise flag
        best_dups.append(choose_best_read(v))
    return best_dups

def choose_best_read(reads):
    best = None
    for r in reads:
        if has_recal_tag(r):
            best = r
            break
        if has_oc_tag(r):
            best = r
    if best is None:
        best = reads[0]
    return best

def has_recal_tag(read):
    return (read.has_tag('BD') or read.has_tag('BI'))
    
def has_oc_tag(read):
    return read.has_tag('OC') 

def get_coordinate(rsum):
    try:
        pos = "{:,}".format(int(rsum.split[3]))
    except ValueError:
        pos = rsum.split[3]
    return "{}:{}".format(rsum.split[2], pos)

def check_header(bamfile):
    header = bamfile.header
    if 'HD' in header:
        if 'SO' in header['HD']:
            if header['HD']['SO'] != 'queryname':
                sys.stderr.write("WARNING: Input does not appear to be sorted"+
                                 " by read name. Header sort order was '{}'\n"
                                  .format(header['HD']['SO']) )
            return
    sys.stderr.write('WARNING: Could not confirm sort order from input header.')

def filter_dups(bam, output, dup_bam=None, ref=None, force=False):
    for fn in (output, dup_bam):
        if fn is not None and os.path.exists(fn) and not force:
            sys.exit("Output file '{}' exists - please delete or choose"
                     .format(fn) + " another name.")
    bamfile = get_bamfile(bam)
    outfile = get_output(output, bamfile) 
    dupfile = None
    if dup_bam is not None:
        dupfile = get_output(dup_bam, bamfile)
    n = 0
    d = 0
    f = 0
    prev_id = None
    prog_string = ''
    cache = []
    check_header(bamfile)
    for read in bamfile.fetch(until_eof=True):
        n += 1
        if prev_id != read.query_name and prev_id is not None:
            valid, dups, ndup = check_dups(cache)
            d += ndup
            f += len(cache) - len(valid)
            for v in valid:
                outfile.write(v)
            if dupfile is not None:
                for dp in dups:
                    dupfile.write(dp)
            cache[:] = []
        cache.append(read)
        prev_id = read.query_name
        if not n % 10000:
            sys.stderr.write("\r" + " " * len(prog_string))
            prog_string = ("{:,} read, {:,} duplicates, {:,} reads filtered "
                           .format(n, d, f)
                          )
            sys.stderr.write("\r" + prog_string)
    valid, dups, ndup = check_dups(cache)
    d += ndup
    f += len(cache) - len(valid)
    for v in valid:
        outfile.write(v)
    if dupfile is not None:
        for dp in dups:
            dupfile.write(dp)
    cache[:] = []
    sys.stderr.write("\nFinished: {:,} alignments read, {:,} reads duplicated," 
                     .format(n, d) + " {:,} reads filtered\n".format(f))
    outfile.close()
    if dupfile is not None:
        dupfile.close()

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    filter_dups(**vars(args))

