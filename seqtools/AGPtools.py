# https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
AGPSEQCODES='WADFGOP'
AGPGAPCODES='NU'
AGPGAPTYPES=['scaffold', 'contig', 'centromere', 'short_arm', 'heterochromatin', 'telomere', 'repeat', 'contamination']
AGPLINKAGES=['yes', 'no']
AGPORIENTATIONS=['+', '-', '?', '0', 'na'] ## By default, components with unknown orientation (?, 0 or na) are treated as if they had + orientation.
AGPEVIDENCE=['na', 'paired-ends', 'align_genus', 'align_xgenus', 'align_trnscrpt', 'within_clone', 'clone_contig', 'map', 'pcr', 'proximity_ligation', 'strobe', 'unspecified']
AGPREORIENT = {'+':'-', '-':'+'}

class AGP_RECORD(object):
    def __init__(self, l):
        '''
        Converts AGP lines into useful object.

        INPUTS:
        l = a list of all the elements of a line from AGP file (i.e. line.strip().split())
        
        '''
        self.agp_record = l
        self.original = {} # for keeping track of orignal state when modifications made

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

    def is_sequence(self):
        return self.component_type() in AGPSEQCODES

    def is_gap(self):
        return self.component_type() in AGPGAPCODES

    def _previously_modified(self, key):
        return key in self.original

    def _original_state(self, key):
        return self.original[key]

    def _add_to_original(self, key, value):
        self.original[key] = value

    def reorient(self):
        assert self.is_sequence()
        if self._previously_modified('orientation'):
            sys.stderr.write("This record's orientation has already been modified. Consider investigating AGP file or code since this should not be possible. Skipping....\n")
        else:
            self._add_to_original('orientation', self.contig_orientation())
            self.agp_record[8] = AGPREORIENT[ self.agp_record[8] ]

    def format_as_line(self):
        return '\t'.join(self.agp_record)

    

def get_ctg_msg(agp_record, verbose=False):
    if verbose:
        scf = agp_record.scaffold()
        ctg = agp_record.contig_id()
        stdsym = agp_record.contig_orientation()
        stdadd = "-" if stdsym == "-" else "+"
        msg ="Building:" + scf + "\n\tAdding:" + ctg + "\n\tStrandSymbol:" + stdsym + "\n\tStrandAdded:" + stdadd + "\n"
        msg += "\t" + " ".join([e for e in agp_record.agp_record]) + "\n"
        sys.stderr.write(msg)


def get_gap_msg(agp_record, verbose=False):
    if verbose:
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
        
def get_scaffold_sequence(scfinfo, contigs, revcomp, verbose=False):
    seq = ''
    for agp_record in scfinfo:
        if agp_record.is_contig():
            get_ctg_msg(agp_record, verbose)
            if agp_record.add_rev_comp():
                seq += revcomp[agp_record.contig_id()]
            else:
                seq += contigs[agp_record.contig_id()]
        elif agp_record.is_gap():
            get_gap_msg(agp_record, verbose)
            if agp_record.gap_length_is_known():
                seq += 'N' * agp_record.gapLength()
            elif agp_record.gap_length_is_unknown():
                seq += 'n' * agp_record.gapLength()
            else:
                error_message("Encountered unexpected gap length status. Only takes N and U gap components for now.")
        else:
            error_message("Encountered unexpected component type. Only takes N, U, and W for now. You can't fire me if I QUIT.")
    return seq

