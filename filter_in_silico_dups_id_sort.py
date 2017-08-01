#!/usr/bin/env python
import sys
import os
from subprocess import PIPE, Popen
import argparse
from collections import defaultdict
import re

oc_re = re.compile(r'OC:Z:(\S+)')
bi_re = re.compile(r'BI:Z:(\S+)')
bd_re = re.compile(r'BD:Z:(\S+)')

class ReadSum(object):
    ''' 
        Class for storing variables for SAM records for convenience
        and readability. Simply stores the binary read, a list of 
        each SAM column as strings, and the bitwise flag as an int.
    '''

    __slots__ = ['rid', 'read', 'split', 'flag', 'is_header']

    def __init__(self, bline):
        '''
            Requires single line of SAM from samtools piped process as 
            bytes
        '''
        self.is_header = False
        self.read = bline
        self.split = bline.decode(sys.stdout.encoding).split()
        if self.split[0][0] == '@': #header
            self.is_header = True
            self.flag = None
            self.rid = None
        else:
            self.flag = int(self.split[1])
            self.rid = self.split[0]
        

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

    required_args = parser.add_argument_group('Required Arguments')
    optional_args = parser.add_argument_group('Optional Arguments')
    required_args.add_argument('bam', metavar='INPUT', 
                               help=''' Input BAM filename''')
    optional_args.add_argument('-o', '--output', metavar='OUT', 
                               help=''' Write cleaned output BAM to this file. 
                                        Default is to write SAM format to 
                                        STDOUT''')
    optional_args.add_argument('-d', '--dup_bam', metavar='DUPS', 
                               help=''' Optional BAM file for duplicated reads.
                               ''')
    optional_args.add_argument('-t', '--threads', metavar='THREADS', default=0,
                               type=int,
                               help=''' Number of additional threads to use 
                                        with samtools sort. ''')
   # optional_args.add_argument('-r', '--ref', metavar='REF', 
   #                            help=''' Reference fasta file. May be required   
   #                                     if using CRAM input.''')
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
        Identify duplicate reads in cache or ReadSummary objects. 

        Returns two lists of binary representations of reads. The first
        is a list of non-duplicated reads, the second is a list of 
        filtered duplicates.
    '''
    non_dups = []
    dups = []
    read1s = []
    read2s = []
    singletons = []
    ndups = 0
    for r in cache:
        if r.flag & 1: #paired
            if r.flag & 64:       #read1
                read1s.append(r)
            elif r.flag & 128:    #read2
                read2s.append(r)
            else:
                sys.exit("ERROR: Malformed BAM - read is paired but " + 
                         "is neither set as first or second in pair for " +
                         "read {}".format(read.decode(sys.stdout.encoding))
                        )
        else:
            singletons.append(r)
    if singletons:
        singles, empty, supp = assign_primary_and_supplementary(singletons)
        ndups += add_non_dups_and_dups(singletons, non_dups, dups)
    r1s, r2s, psupp = assign_primary_and_supplementary(read1s, read2s)
    if psupp:
        ndups += add_non_dups_and_dups(psupp, non_dups, dups)
    if len(r1s) > 1:
        dups.extend(list(x.read for x in r1s))
        ndups += 1
    if len(r2s) > 1:
        dups.extend(list(x.read for x in r2s))
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
        return (r1s[0].read, r2s[0].read)
    pairs = []
    for r1 in r1s:
        for r2 in r2s:
            if (r1.flag & 4 and not r2.flag & 4 or  # unmapped and mapped 
                r2.flag & 4 and not r1.flag & 4):   # - expect same coordinate,
                #coords may differ if mapped read underwent indel realignment
                #and unmapped read originated from before indel realignment
                if r1.split[3] == r2.split[3] and r1.split[2] == r2.split[2]:
                    pairs.append((r1, r2))
            elif not r1.flag & 4 and not r2.flag & 4: #both mapped 
                #check mate coords match
                if (r1.split[3] == r2.split[7] and 
                    (r1.split[2] == r2.split[6] or r2.split[6] == '=')):
                    pairs.append((r1, r2))
            else:#both unmapped
                pairs.append((r1, r2))
    best_pair = None
    for p in pairs:
        if has_recal_tag(p[0]) and has_recal_tag(p[1]):
            best_pair = (p[0].read, p[1].read)
            break
        elif has_oc_tag(p[0]) or has_oc_tag(p[1]):
            best_pair = (p[0].read, p[1].read)
    if best_pair is None and pairs:
        best_pair = (pairs[0][0].read, paris[0][1].read)
    return best_pair

def add_non_dups_and_dups(reads, non_dups, dups):
    sdups,snon_dups = get_dup_reads(reads)
    non_dups.extend(snon_dups)
    if sdups:
        non_dups.extend(choose_best_dups(sdups))
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
                if reads[i].flag & 256 or reads[i].flag & 2048: #not primary/supplementary 
                    if isplit[2] == jsplit[2] and isplit[3] == jsplit[3]:
                        is_dup = True
                        dup_sets[reads[i].flag].append(reads[i])
                else:       
                    is_dup = True
                    dup_sets[reads[i].flag].append(reads[i])
            if is_dup:
                break
        if not is_dup:
            non_dups.append(iread)
    return dup_sets, non_dups

def choose_best_dup(dups):
    best_dups = []
    for k,v in dups.items(): #key is bitwise flag
        best_dups.append(choose_best_read(read))

def choose_best_read(reads):
    best = None
    for r in reads:
        if has_recal_tag(r):
            best = r.read
            break
        if has_oc_tag(r):
            best = r.read
    if best is None:
        best = reads[0].read
    return best

def has_recal_tag(read):
    for tag in read.split[11:]:
        if bi_re.match(tag) or bd_re.match(tag):
            return True
    
def has_oc_tag(read):
    for tag in read.split[11:]:
        if oc_re.match(tag):
            return True

def get_coordinate(rsum):
    try:
        pos = "{:,}".format(int(rsum.split[3]))
    except ValueError:
        pos = rsum.split[3]
    return "{}:{}".format(rsum.split[2], pos)

def filter_dups(bam, output=None, dup_bam=None, ref=None, threads=0):
    for fn in (output, dup_bam):
        if fn is not None and os.path.exists(fn):
            sys.exit("Output file '{}' exists - please delete or choose"
                     .format(fn) + " another name.")
    if output is None:
        output = sys.stdout
    outfile = open(output, 'wb')
    dupfile = None
    dupwrite = None
    dash = '-' #subprocess doesn't seem to like '-' passed as a string
    samwrite = Popen(["samtools", 'view', '-Sbh', dash], stdout=outfile, 
                     stdin=PIPE, bufsize=-1)
    if dup_bam is not None:
        dupfile = open(dup_bam, 'wb')
        dupwrite = Popen(["samtools", 'view', '-Sbh', dash], stdout=dupfile,
                         stdin=PIPE, bufsize=-1)
    n = 0
    d = 0
    f = 0
    prev_id = None
    prog_string = ''
    cache = []
    sys.stderr.write("Sorting input file by read ID... This may take some " + 
                     "time...\n")
    samsort = Popen(["samtools", 'sort', '-@', str(threads), '-n', '-O', 'SAM', 
                     bam], 
                    stdout=PIPE, bufsize=-1)
    for line in samsort.stdout:
        rsum = ReadSum(line)
        if rsum.is_header:
            samwrite.stdin.write(rsum.read)
            if dupwrite is not None:
                dupwrite.stdin.write(rsum.read)
            continue
        n += 1
        if prev_id != rsum.rid and prev_id is not None:
            valid, dups, ndup = check_dups(cache)
            d += ndup
            f += len(cache) - len(valid)
            for v in valid:
                samwrite.stdin.write(v)
            if dupwrite is not None:
                for dp in dups:
                    dupwrite.stdin.write(dp)
            cache[:] = []
        cache.append(rsum)
        prev_id = rsum.rid
        if not n % 10000:
            sys.stderr.write("\r" + " " * len(prog_string))
            prog_string = ("{:,} read, {:,} duplicates, {:,} reads filtered "
                           .format(n, d, f)
                          )
            sys.stderr.write("\r" + prog_string)
    sys.stderr.write("\nFinished: {:,} alignments read, {:,} reads duplicated," 
                     .format(n, d) + " {:,} reads filtered\n".format(f))
    if output is not None:
        outfile.close()
    if dupfile is not None:
        dupfile.close()

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    filter_dups(**vars(args))

