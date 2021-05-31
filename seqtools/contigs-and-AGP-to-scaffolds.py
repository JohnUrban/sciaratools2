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
Takes in Contig Fasta and AGP.

OUPUT:
Outputs scaffolds fasta file.

ASSUMPTIONS:
1. Names (record.id) of contig sequences in FASTA match AGP names in column 6.
2. Lines in AGP are ordered for each scaffold - i.e. not all mixed up and out of order for w/e reason that would be.


OTHER:
AGP coordinates are 1-based inclusive.

AGP entry for contig:
name start end part_number component_type name_id contig_start contig_end orientation
component_type = W
orientation = +


AGP entry for gap:
name start end part_number component_type gapLen gapType linkage evidence
component_type = N
gapType = scaffold
linkage = yes
evidence = map


Example of scaffold with 2 contigs separated by one gap from map-evidence:
trial_scaffold_1        1       9       1       W       trial_scaffold_1_1      1       9       +
trial_scaffold_1        10      33      2       N       24      scaffold        yes     map
trial_scaffold_1        34      39      3       W       trial_scaffold_1_3      1       6       +


Updates:
5/28/2021   It now accepts AGP files with U components, in addition to N and W.
            It writes N gaps of known length with uppercase NNNN and U gaps of unknown length with lowercase nnnn.
            It is now more orientation-aware, and reverse-complements contigs when needed before adding to a growing scaffold.
            Previously, my only use cases did not require this.
            In total, this now works with AGP files provided by Phase Genomics for Hi-C scaffolds.
            It also (already) works with AGP files produced from any set of scaffolds (e.g. bioNano) produced with scf-N-to-AGP.py.
            Note: scf-N-to-AGP.py is now also updated to allow creation of U gap entries.
            This script also can now handle comment/header lines starting with "#" (as in Phase Genomics AGP files).

''', formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("fasta")
parser.add_argument("-a", "--agp", type=str, required=True, help='''Path to AGP.''')
parser.add_argument("-v", "--verbose", action='store_true', default=False, help='''Messages to stderr about progress.''')
args = parser.parse_args()


## DEFINITIONS
## Not "from AGPtools import *"


#########################################################
## Read in FASTA into dictionary
contigs = {}
revcomp = {}
for record in SeqIO.parse(args.fasta, 'fasta'):
    contigs[str(record.id)] = str(record.seq)
    revcomp[str(record.id)] = str(record.seq.reverse_complement())
    #ans = contigs[str(record.id)] == revcomp[str(record.id)]


## Read in AGP:
AGP = open(args.agp).readlines()

## Go through AGP and create defaultdict that maps scaffold names to agp lines
SCF = defaultdict(list)
ORDER = []
for line in AGP:
    if line.startswith('#'):
        continue
    line = line.strip().split()
    SCF[line[0]].append( AGP_RECORD(line) )

    if len(ORDER) == 0 or ORDER[-1] != line[0]:
        ORDER.append( line[0] )

## Go through scaffold information and put together
for scaffold in ORDER:
    # Return FASTA name
    print(">"+scaffold)

    # Return FASTA sequence
    print(get_scaffold_sequence(SCF[scaffold], contigs, revcomp, args.verbose))


        
        



# END


