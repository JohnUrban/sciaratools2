#!/usr/bin/env python3

import argparse, sys
import re
from Bio import SeqIO
from collections import defaultdict

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
class AGP_RECORD(object):
    def __init__(self, l):
        '''
        Converts AGP lines into useful object.

        INPUTS:
        l = a list of all the elements of a line from AGP file (i.e. line.strip().split())
        
        '''
        self.agp_record = l

    def scaffold(self):
        return self.agp_record[0]

    def scaffold_start(self, astype=int):
        return astype(self.agp_record[1])

    def scaffold_end(self, astype=int):
        return astype(self.agp_record[2])

    def part_number(self, astype=int):
        return astype(self.agp_record[3])

    def component_type(self, astype=str):
        return astype(self.agp_record[4])

    def contig_id(self, astype=str):
        return astype(self.agp_record[5])

    def gapLength(self, astype=int):
        return astype(self.agp_record[5])

    def contig_start(self, astype=int):
        return astype(self.agp_record[6])

    def contig_end(self, astype=int):
        return astype(self.agp_record[7])

    def contig_orientation(self, astype=str):
        return astype(self.agp_record[8])

    def gapType(self, astype=str):
        return astype(self.agp_record[6])

    def linkage(self, astype=str):
        return astype(self.agp_record[7])

    def evidence(self, astype=str):
        return astype(self.agp_record[8])

    def is_gap(self):
        return self.component_type() == 'N' or self.component_type() == 'U' 

    def is_contig(self):
        return self.component_type() == 'W'

    def add_rev_comp(self):
        return self.contig_orientation() == "-"

    def gap_length_is_known(self):
        return self.component_type() == 'N'
    
    def gap_length_is_unknown(self):
        return self.component_type() == 'U' 

def get_ctg_msg(args, agp_record):
    if args.verbose:
        scf = agp_record.scaffold()
        ctg = agp_record.contig_id()
        stdsym = agp_record.contig_orientation()
        stdadd = "-" if stdsym == "-" else "+"
        msg ="Building:" + scf + "\n\tAdding:" + ctg + "\n\tStrandSymbol:" + stdsym + "\n\tStrandAdded:" + stdadd + "\n"
        msg += "\t" + " ".join([e for e in agp_record.agp_record]) + "\n"
        sys.stderr.write(msg)


def get_gap_msg(args, agp_record):
    if args.verbose:
        scf = agp_record.scaffold()
        gaplenstatus = "Known" if agp_record.gap_length_is_known() else "Unknown"
        gapchar = "N" if agp_record.gap_length_is_known() else "n"
        msg ="Building:" + scf + "\n\tAdding:Gap\n\tGapLength:" + str(agp_record.gapLength()) + "\n"
        msg += "\tGapLengthStatus:" + gaplenstatus + "\n\tGapCharacter:" + gapchar + "\n"
        msg += "\t" + " ".join([e for e in agp_record.agp_record]) + "\n"
        sys.stderr.write(msg)

def error_message(msg=""):
    sys.stderr.write("ERROR: " + msg + "\n")
    sys.stderr.write("You can't fire me if I QUIT.\n")
    quit()
        
def get_scaffold_sequence(scfinfo, contigs, revcomp):
    seq = ''
    for agp_record in scfinfo:
        if agp_record.is_contig():
            get_ctg_msg(args, agp_record)
            if agp_record.add_rev_comp():
                seq += revcomp[agp_record.contig_id()]
            else:
                seq += contigs[agp_record.contig_id()]
        elif agp_record.is_gap():
            get_gap_msg(args, agp_record)
            if agp_record.gap_length_is_known():
                seq += 'N' * agp_record.gapLength()
            elif agp_record.gap_length_is_unknown():
                seq += 'n' * agp_record.gapLength()
            else:
                error_message("Encountered unexpected gap length status. Only takes N and U gap components for now.")
        else:
            error_message("Encountered unexpected component type. Only takes N, U, and W for now. You can't fire me if I QUIT.")
    return seq




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
    print(get_scaffold_sequence(SCF[scaffold], contigs, revcomp))


        
        



# END


