#!/usr/bin/env python3

import argparse, sys
import re
from Bio import SeqIO

# Parse command-line arguments
parser = argparse.ArgumentParser(description='''
Description: Takes in Fasta.
Checks beginning and end of each sequence to ensure it does not start with N (or w/e specified).
''')
parser.add_argument("fasta")
parser.add_argument("-G", "--gapchar", type=str, default='N', help='''The single letter gap character. Default =  N.''')
parser.add_argument("-i", "--ignorecase", action='store_true', default=False, help='''When given a single gap character, this tells it to look for upper and lower case.''')
args = parser.parse_args()


# Open FASTA, search for masked regions, print in BED3 format
if args.ignorecase:
    N = '['+ args.gapchar.upper() + args.gapchar.lower() + ']+'
else:
    N = args.gapchar + '+' ## 'N+'
handle = sys.stdin if args.fasta in ('-','stdin') else open(args.fasta)
pattern = N

## HEADER
contents = ['#', 'name', 'length', 'leading_N', 'trailing_N', 'first_N_pos', 'last_N_pos', 'num_Intervals', 'length_N', 'percent_N', 'length_ACGT']
print('\t'.join(contents))

## RECORDS
for record in SeqIO.parse(handle, "fasta"):
    length = len(record.seq)
    allmatch = [e for e in re.finditer(pattern, str(record.seq))]
    nIntervals = len(allmatch)
    nN = 0
    if nIntervals > 0:
        first = allmatch[0].start()
        last = allmatch[-1].end()
        start = 'yes' if first == 0 else 'no'
        end = 'yes'if last == length else 'no'
        for e in allmatch:
            nN += e.end()-e.start()
        #contents = [record.id, length, start, end, first, last, nN, pctN]
    else: 
        start = 'no'
        end = 'no'
        first = 'na'
        last = 'na'
        #contents = [record.id, length, 'no','no', 'na', 'na', nN, pctN]
    pctN = round(100.0*nN/length, 3)
    nACGT = length - nN
    contents = [record.id, length, start, end, first, last, nIntervals, nN, pctN, nACGT]
    print('\t'.join([str(e) for e in contents]))

    
if args.fasta not in ('-','stdin'):
    handle.close()
