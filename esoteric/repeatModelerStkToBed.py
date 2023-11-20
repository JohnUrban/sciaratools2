#!/usr/bin/env python3

import sys, argparse
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser(description="""

Created 2023-11-20.
	Latest update on:	2023-11-20 
	Other updates on:	2023-11-20

Take in a STK file from RepeatModeler.

Return BED file of seed locations for each repeat family.


NOTES:
    - STK is in Stockholm format used by Dfam, Rfam, pFam
        - https://en.wikipedia.org/wiki/Stockholm_format
        - https://github.com/Dfam-consortium/RepeatModeler/tree/master
    - From examples seen, and from exploring my own STK outputs from RM, the seq coords given are 1-based.
    - Coords are given stranded... so if on neg strand, they are backwards.



    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--input', '-i',
                   type=str, required=True,
                   help='''Path to input STK file''')
parser.add_argument('--outpre', '-o',
                   type=str, default='stdout',
                   help='''Prefix for output files. Default: stdout.''')

args = parser.parse_args()

########################################################################################################
## FXNS

def stkToBed(stk_fh, outpre=None):
    assert outpre is not None
    bedout = sys.stdout if outpre == 'stdout' else open(outpre+'.bed','w')
    bedout.write(('\t').join(['#seq', 'start', 'end', 'ID(famiy name);TP(type);CC(reflen)','SQ','.','alignment'])+'\n')
        
    with open(stk_fh) as stk:
        family = None
        info = {}
        for line in stk:
            line = line.strip()
            headerline = False
            if line.startswith('#'):
                headerline = True
                line = line.split()
                info[line[1]] = (' ').join(line[2:])

            emptyline = line == ''
            endoffile = line == '//'
            
            if not headerline and not emptyline and not endoffile:
                
                ##BED-like
                line = line.split()
                bed3 = line[0].split(':')
                seq = bed3[0]
                coords = [int(e) for e in bed3[1].split('-')]
                if coords[0] < coords[1]:
                    start = str(coords[0]-1)
                    end = coords[1]
                else:
                    start = str(coords[1]-1)
                    end = coords[0]
                aln = line[1]
                reflen = info['CC'].split()[0].split('=')[1]+'bp' 
                bedlist = [seq, start, end, info['ID']+';'+info['TP']+';'+reflen, info['SQ'],'.',aln]
                           
          
                bedline = ('\t').join([str(e) for e in bedlist])
                bedout.write( bedline + '\n' )
                      
                    
    if outpre != 'stdout':
        bedout.close()
    return None




try:                
    stkToBed(args.input, outpre=args.outpre)          
except (BrokenPipeError, IOError):
    sys.stderr.close()

########################################################################################################
## BODY
