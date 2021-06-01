from Bio import SeqIO
import sys
from collections import defaultdict

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

    def renumber(self, newnumber):
        self.agp_record[3] = str( newnumber )

    def update_scaffold_start(self, newstart):
        self.agp_record[1] = str( newstart )

    def update_scaffold_end(self, newend):
        self.agp_record[2] = str( newend )

    def format_as_line(self):
        return '\t'.join(self.agp_record)


class AGP_FILE(object):
    def __init__(self, AGP):
        ''' AGP     =   string, agp file path.'''
        # store filename
        self.filename = AGP

        # Header storage variable 
        self.header = []

        # Scaffold-centric storage dict for AGP records
        self.SCF = defaultdict(list)

        # List variable to store original scaffold ordering in AGP file (lost in dict)
        self.ORDER = []

        # pre-process agp lines
        AGP = open(self.filename).readlines() ## not "self.AGP" on purpose.

        ## Go through AGP and:
        ##  - map scaffold names to list of all agp lines corresponding to that scaffold in order of appearance.
        ##  - keep track of scaffold order of appearance.
        for line in AGP:
            if line.startswith('#'):
                self.header.append(line)
                continue
            line = line.strip().split()
            self.SCF[line[0]].append( AGP_RECORD(line) )

            if len(self.ORDER) == 0 or self.ORDER[-1] != line[0]:
                self.ORDER.append( line[0] )

    def get_ordered_scf_list(self):
        return self.ORDER

    def get_agp_records_for_scf(self, scf_name):
        return self.SCF[scf_name]

    def reverse_order_of_records_in_scf(self, scf_name):
        ''' This will reverse the order of appearance, but will not adjust the +/- signs for contigs for full reverse complement.'''
        self.SCF[scf_name] = list( reversed( self.SCF[scf_name] ) )
        self.update_part_numbers_in_scf(scf_name)
        self.update_coordinates_in_scf(scf_name)

    def update_part_numbers_in_scf(self, scf_name):
        '''
        Will go through and re-number the parts in order they appear.
        Use case = after reversing order of parts in scf.
        '''
        newnumber = 0
        for agp_record in self.get_agp_records_for_scf(scf_name):
            newnumber += 1
            agp_record.renumber( newnumber )

    def update_coordinates_in_scf(self, scf_name):
        '''
        Will go through and re-write scf start/end given the order components appear.
        Use case = after reversing order of parts in scf.
        '''
        newstart = None
        newend = None
        for agp_record in self.get_agp_records_for_scf(scf_name):
            if agp_record.is_sequence():
                a = agp_record.contig_start()
                b = agp_record.contig_end()
                component_length = b - a + 1
            elif agp_record.is_gap():
                component_length = agp_record.gapLength()
            else:
                error_message("Malformed AGP record encountered....")

            if newstart is None and newend is None:
                newstart = 1
                newend = component_length
            else:
                newstart = newend + 1
                newend = newend + component_length

            agp_record.update_scaffold_start( newstart )
            agp_record.update_scaffold_end( newend )


    def reorient_all_records_in_scf(self, scf_name):
        for agp_record in self.get_agp_records_for_scf(scf_name):
            if agp_record.is_sequence():
                agp_record.reorient()

    def reverse_complement_full_scf(self, scf_name):
        ''' This will reverse the order of appearance of records AND adjust the +/- signs for contigs'''
        self.reverse_order_of_records_in_scf(scf_name)
        self.reorient_all_records_in_scf(scf_name)

    def reverse_complement_full_scf_for_all_named(self, names, verbose = False):
        ''' This will reverse the order of appearance of records AND adjust the +/- signs for contigs'''
        success = [] ## For keeping tabs on names found
        for scf_name in self.get_ordered_scf_list():
            if scf_name in names:
                self.reverse_complement_full_scf( scf_name )
                success.append( scf_name )
        if verbose:
            self._warn_names(names, success)
        
        
    def reorient_named_contigs_in_scf(self, scf_name, names):
        success = [] ## For keeping tabs on names found
        for agp_record in self.get_agp_records_for_scf(scf_name):
            if agp_record.is_sequence():
                if agp_record.contig_id() in names:
                    agp_record.reorient()
                    success.append( agp_record.contig_id() )
        return success

    def reorient_named_contigs_in_AGP(self, names, verbose = False):
        ''' In any Scaffold.'''
        success = [] ## For keeping tabs on names found (need to do list(set(l)) later)
        for scf_name in self.ORDER:
            success += self.reorient_named_contigs_in_scf(scf_name, names)
        if verbose:
            self._warn_names(names, success)

    def write_out(self, fhconn, include_header=True):
        ''' Assumes fhconn is stdout or an open file connection)'''
        ## Write header if present, and desired.
        if include_header and len(self.header) > 0:
            for hline in self.header:
                fhconn.write(hline)
        ## Write out AGP record lines
        for scf in self.get_ordered_scf_list():
            for agp_record in self.get_agp_records_for_scf( scf ):
                fhconn.write( agp_record.format_as_line() + '\n' )
                

    def _warn_names(self, names, success):
        names = set(names)
        success = set(success)
        failed = list( names.difference(success) )
        n = len(names)
        f = len(failed)
        warn_message("There were " + str(f) + " of " + str(n) + ' names not found :')
        for name in names:
            warn_message("\t- " + name + " not found....")



def reorient(AGP, names, target, verbose=False):
    if target == 's':
        AGP.reverse_complement_full_scf_for_all_named(names, verbose)
    elif target == 'c':
        AGP.reorient_named_contigs_in_AGP(names, verbose)
    else:
        error_message("Encountered unanticipated target code....")
        
def write_out_AGP_file_from_object(AGP, fh, include_header=True):
    ## Open output connection
    outhandle = open_fh_connection( fh, 'w' )

    # Write out modified AGP
    AGP.write_out( outhandle,
                   include_header = include_header )

    # Close output connection
    close_fh_connection( fh, outhandle )


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

def error_message(msg="", fireself=True):
    sys.stderr.write("ERROR: " + msg + "\n")
    if fireself:
        sys.stderr.write("You can't fire me if I QUIT.\n")
        quit()

def warn_message(msg="", fireself=False):
    sys.stderr.write("WARNING: " + msg + "\n")
    if fireself:
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


def build_contig_dicts(fasta):
    ## Read in FASTA into dictionary
    contigs = {}
    revcomp = {}
    for record in SeqIO.parse(fasta, 'fasta'):
        contigs[str(record.id)] = str(record.seq)
        revcomp[str(record.id)] = str(record.seq.reverse_complement())
        #ans = contigs[str(record.id)] == revcomp[str(record.id)]
    return contigs, revcomp

def open_fh_connection(fh, mode='w'):
    if mode == 'w':
        return sys.stdout if fh in ('-', 'stdout', sys.stdout, None) else open(fh, mode)
    elif mode == 'r':
        return sys.stdin if fh in ('-', 'stdin', sys.stdin, None) else open(fh, mode)
    else:
        error_message("Unanticipated file mode encountered....")

def close_fh_connection(fh, fhconn):
    if fh not in ('-', 'stdin', 'stdout', 'stderr', sys.stdout, sys.stdout, sys.stderr, None):
        fhconn.close()
        
def build_and_print_scaffolds_from_contigs_and_agp(contig_fasta, agp, outfname=sys.stdout, verbose=False):
    ## Read in and build contig dicts
    contigs, revcomp = build_contig_dicts(contig_fasta)

    ## Read in and build AGP object
    AGP = AGP_FILE(agp)

    ## Open connection to output
    out = open_fh_connection(outfname)

    ## Go through scaffold information and put together
    for scaffold in AGP.get_ordered_scf_list():
        # Return FASTA name
        out.write('>' + scaffold + '\n')

        # Return FASTA sequence
        out.write(
            get_scaffold_sequence(
                AGP.get_agp_records_for_scf(scaffold),
                contigs,
                revcomp,
                verbose)
            + '\n')

    ## Close file connection if not written to stdout
    close_fh_connection(outfname, out)




def is_target_for_modification(agp_record, names):
    in_names = (agp_record.scaffold() in names or agp_record.contig_id() in names)
    return in_names

def get_line_list(fh):
    handle = open_fh_connection(fh, 'r')
    lines = [line.strip() for line in handle]
    close_fh_connection(fh, handle)
    return lines




