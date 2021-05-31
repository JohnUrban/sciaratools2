#!/usr/bin/env python3

import argparse, sys
import re
from Bio import SeqIO
from collections import defaultdict
from AGPtools import *

# Parse command-line arguments
parser = argparse.ArgumentParser(description='''
Description:

INPUT:
Takes in AGP and list of names (found in either col1 or col6 of AGP) for modifying.

OUPUT:
AGP file modified according to instructions.
- E.g. --reorient will change + to - or - to +.

NOTE:
In current version,
- it will look for given names in both column 1 and 6 (scaffnames and contig names respectively).
- it will make modification if name is in either one.
- this is fine for now since there should be no multi-contig scaffold where the scaffold has the same name as one contig.
    - single-contig scaffolds with the same name as the single contig will follow expected behavior as well
- If name is a scaffold name, then every contig on the scaffold will be re-oriented (or otherwise modified).
- If name is a contig name, only that contig will be re-oriented (or otherwise modified).

OTHER:
AGP coordinates are 1-based inclusive.

AGP entry for contig:
name start end part_number component_type name_id contig_start contig_end orientation

Typical examples:
- component_type = W
- orientation = +


AGP entry for gap:
name start end part_number component_type gapLen gapType linkage evidence

Typical examples:
- component_type = N or U
- gapType = scaffold
- linkage = yes
- evidence = map or paired-ends


Example of scaffold with 2 contigs separated by one gap from map-evidence:
trial_scaffold_1        1       9       1       W       trial_scaffold_1_1      1       9       +
trial_scaffold_1        10      33      2       N       24      scaffold        yes     map
trial_scaffold_1        34      39      3       W       trial_scaffold_1_3      1       6       +


''', formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-a", "--agp", type=str, required=True, help='''Path to AGP.''')
parser.add_argument("-n", "--names", type=str, default=False, help='''Path single-column file with one name per line. Can be from stdin (use: -, stdin).''')
parser.add_argument("-c", "--namescmd", type=str, default=False, help='''Comma-separated list from commandline.''')
parser.add_argument("-R", "--reorient", action='store_true', default=False, help='''Re-orient sequences. If +, then -. If -, then +.''')
parser.add_argument("-o", "--outputname", type=str, default=None, help='''AGP output filename. Default is to stdout.''')

args = parser.parse_args()

## Definitions

def is_target_for_modification(agp_record, names):
    in_names = (agp_record.scaffold() in names or agp_record.contig_id() in names)
    return in_names


## Get names:
names = []
if args.names:
    #OPEN CONNECTION TO FILE/STDIN 
    handle = sys.stdin if args.names in ('-','stdin') else open(args.names)
    names += [line.strip() for line in handle]
    if args.names not in ('-','stdin'):
        handle.close()
if args.namescmd:
    names += [e.strip() for e in args.namescmd.strip().split(',')]


## Read in AGP:
AGP = [AGP_RECORD(line.strip().split()) for line in open(args.agp).readlines() if not line.startswith('#')]

## Open output connection
outhandle = sys.stdout if args.outputname is None else open(args.outputname, 'w')


## Go through
for agp_record in AGP:
    # Make modifications according to options flagged IF record is a target for modifications.
    if is_target_for_modification(agp_record, names):
    
        # 1. Re-orientation
        if args.reorient and agp_record.is_sequence():
            agp_record.reorient()

        # 2. Something else... etc.
    
    # Return line / Write out
    outhandle.write( agp_record.format_as_line() + '\n' )
        
# Close output connection
if args.outputname is not None:
    outhandle.close()
