#!/usr/bin/env python3

import sys, argparse
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser(description="""

Created 2023-11-17.
	Latest update on:	2023-11-17 
	Other updates on:	2023-11-17

Take in a DAT file from TRF.

Return other stuff.

- BED file of TR locations and info (copy number, etc)




TODO/ADD MAYBE.
- FASTA file of TR units with location info (gives array locations).
    - Can be made from the BED-like output outside this program.
- FASTA file of TR arrays with location info (and unit info)
    - Can be made from the BED-like output outside this program.
- statistics for each sequence
    - how many TRs
    - how many unique TRs
- statistics for each unique TR unit
    - how many times are arrays found with this unit
    - how many total copies of that unit across all arrays (sum)
    - min, max, mean, median, etc for arrays with that unit
    - Can use the BED file with awk and other already-made programs to get all this.
- are any shorter units found within longer units....
    - for each unit,
    - (or just use BLAST).


    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--input', '-i',
                   type=str, required=True,
                   help='''Path to input DAT file''')
parser.add_argument('--outpre', '-o',
                   type=str, default='stdout',
                   help='''Prefix for output files. Default: stdout.''')

args = parser.parse_args()

########################################################################################################
## FXNS

def datToDictOrBedOrBoth(dat_fh, build_dict=False, outpre=False):
    if build_dict:
        d = defaultdict(list)
    if outpre:
        bedout = sys.stdout if outpre == 'stdout' else open(outpre+'.bed','w')
        bedout.write(('\t').join(['#seq', 'start', 'end', 'period_size','copy_number','consensus_size','percent_matches','percent_indels','score','A','C','G','T','Entropy','Unit','Array'])+'\n')
        
    with open(dat_fh) as dat:
        seq = None
        for line in dat:
            line = line.strip()
            headerline = False
            for word in ['Tandem', 'Program', 'Boston', 'Version','Parameter', 'Sequence']:
                if line.startswith(word):
                    headerline = True
            emptyline = line == ''
            if line.startswith('Sequence'):
                seq = line.split()[1]
                if build_dict:
                    d['seqs'].append( seq )
            if not headerline and not emptyline and seq is not None:
                if build_dict:
                    d[seq].append( line.split() )
                if outpre:
                    ##BED-like
                    line = line.split()
                    bedlist = [seq, str(int(line[0])-1), line[1]] + line[2:]
                               
                    bedline = ('\t').join(bedlist)
                    bedout.write( bedline + '\n' )
                      
                    
    if outpre:
        bedout.close()
    if build_dict:
        return d
    return None




try:                
    datToDictOrBedOrBoth(args.input, build_dict=False, outpre=args.outpre)          
except (BrokenPipeError, IOError):
    sys.stderr.close()
