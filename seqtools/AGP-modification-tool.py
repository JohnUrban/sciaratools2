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
In current version of re-orientation,
- it will only look for given names in column 1 or (XOR) column 6 (scaffnames or contig names respectively).
- If scaffold names are chosen, it will:
    - reverse the order of contig appearance,
    - re-number the scaffold parts s.t. in an n-part scaffold, the formerly nth part is now 1, and the formerly-first part is now n.
    - adjust the coordinates in columns 2 and 3.
    - re-orient strand (+ > -; - > +).
    - NOTE: in downstream applications, compared to upstream ones, re-numbering scaffold parts means it will break relationships relying on identifying a ctg by a scaffold and its part number.
        - I can make it optional to not re-number in the future, but personally have no use cases for that.
- If contig names are chosen, each contig named will be re-oriented in place on the scaffold. 

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
parser.add_argument("-T", "--target", type=str, default="s", required=True, help='''Provide "s" or "c" telling the script whether to focus on scaffolds or contigs only.''')
parser.add_argument("-v", "--verbose", action='store_true', default=False, help='''Say stuff.''')
parser.add_argument("-X", "--exclude_header", action='store_true', default=False, help='''Exclude header from output.''')

args = parser.parse_args()

## Get names:
names = []
names += get_line_list(args.names) if args.names else []
names += [e.strip() for e in args.namescmd.strip().split(',')] if args.namescmd else []


## Read in and build AGP object
AGP = AGP_FILE(args.agp)


# Process AGP according to mod instructions
# 1. Re-orientation
if args.reorient:
    reorient(AGP, names, args.target, args.verbose)
    

# Write out modified AGP
write_out_AGP_file_from_object(AGP = AGP,
                               fh = args.outputname,
                               include_header = not args.exclude_header)

