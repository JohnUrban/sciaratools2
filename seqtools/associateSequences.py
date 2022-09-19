#!/usr/bin/env python3

import argparse, sys, os, shutil, pybedtools, gzip, io
from collections import defaultdict, Counter
import numpy as np

## SCIARA TOOLS IMPORTS
SCITOOLS_BIN_PATH = os.path.dirname(__file__)
if not SCITOOLS_BIN_PATH in sys.path:
    sys.path.append(SCITOOLS_BIN_PATH)
from utilPAF2QueryCentricBEDlike import convertPAFtoBEDlike
from partitionOverlappingBEDintervals import SliceLine, slice3, run_sliceToBed, sortbed, create_sorted_slice_file_from_bed


##############################################################################
''' FUNCTIONS '''
##############################################################################

def parse_args():

    
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter, description='''
    Upstream work:
        Map 1 set of sequences to another (e.g. with Minimap2) and store alignments as PAF. Suggested approach: minimap2 --secondary=no -x asm20 ${primaryseqs} ${assocseqs}

    Description:
        Takes in PAF, and optionally a gap BED file or faCount TSV file, and optional BEDPE.
            gap BED = BED file with intervals for gaps on each associated sequence. Determines N count by summing intervals for each.
            faCount file = output from UCSC/Kent Tools faCount ; assumes N count for each chromosome is in column 7. Can also specify a column if using another format (e.g. seqname Ncount).
        Outputs most likely primary sequence for each associated sequence.
            For PAF, it reports primary with highest percent of assoc mapped, and does this for multiple MAPQ filtering.
            For BEDPE, it reports primary with most links and primary with highest obs:exp ratio of links.
        
    Inputs:
        - PAF
        - optional BED or TSV for N-counts
        - optional BEDPE with links between primary and assoc contigs
        
        
    ''')

    ## Analysis of contig alignments to chromosome scaffolds
    parser.add_argument('-p', '--paf', type=str, required=True,
                        help='''PAF alignment file of assoc seqs mapped to primary seqs.''')


    ncounta = parser.add_mutually_exclusive_group(required=True)
    ncounta.add_argument("-b", "--bed",         type=str, default=False,
                         help='''Optional BED file with N-and/or-n-gap intervals for associated sequences.''')
    ncounta.add_argument("-c", "--facount",     type=str, default=False,
                         help='''Optional Tab-sep text file with seq names and N-counts for associated seqs. Assumes faCount output, seq name in column 1, seq lengths in column 2, and N count in column 7. Use --lengthcol and --ncountcol to change.''')
    ncounta.add_argument("-t", "--seqtable",    type=str, default=False,
                         help='''Same as --facount, but makes no assumptions on length and N count columns. Both -N and -L need to be specified. Script will fail if they are not. Ideal when using any other table to prevent analyzing incorrect columns.''')

    parser.add_argument("-G", "--lengthsfile",  type=str, default=False,
                         help='''Two-column, tab-separated values file with seqnames and lengths for associated seqs. Same as what is called a genome file with bedtools. If BED used rather than faCount output, then this is required. If this given with facount, this supercedes lengths learned from facount.''')
    parser.add_argument("-S", "--seqnamecol",   type=int, default=None,
                        help='''1-based Column index of --facount/--seqtable argument/file to find sequence name.''')
    parser.add_argument("-N", "--ncountcol",    type=int, default=None,
                        help='''1-based Column index of --facount/--seqtable argument/file to find N-counts in.''')
    parser.add_argument("-L", "--seqlengthcol",    type=int, default=None,
                        help='''1-based Column index of --facount/--seqtable argument/file to find sequence lengths in.''')

    parser.add_argument("--min_pct_length_in_matches", type=float, default=0.1,
                        help='''.''')
    parser.add_argument("--min_pct_length_in_alignments", type=float, default=0.2,
                        help='''.''')
    parser.add_argument("--min_prob_given_matches", type=float, default=0.7,
                        help='''.''')
    parser.add_argument("--min_prob_given_alignments", type=float, default=0.7,
                        help='''.''')
    parser.add_argument("--min_mode_of_assignment", type=float, default=0.7,
                        help='''Percent of all tests assignment must pass. Ntests = 4 * num_mapq_vals''')


    parser.add_argument("-m", "--pafmapq", type=str, default='0,10',
                        help='''MAPQ cutoffs to analyze PAF. Provide comma-separated integers. Default = 0,10''')


##    ## Hi-C Linkage analysis
##    parser.add_argument("-B", "--bedpe", type=str, default=False, help='''Optional bedpe file (bamToBed output style).''')
##    ncountp = parser.add_mutually_exclusive_group()
##    ncountp.add_argument("-bp", "--bedp",       type=str, default=False,
##                         help='''Only used with --bedpe option. Optional BED file with N|n-gap intervals for primary sequences. ''')
##    ncountp.add_argument("-cp", "--facountp",   type=str, default=False,
##                         help='''Only used with --bedpe option. Optional Tab-sep text file with seq names and N-counts for primary seqs. Only used with --bedpe option. Assumes faCount output, and N count in column 7. Use --ncol to change.''')
##    ncounta.add_argument("-tp", "--seqtablep",    type=str, default=False,
##                         help='''Only used with --bedpe option. Same as --facountp, but makes no assumptions on length and N count columns. Both -Np and -Lp need to be specified. Script will fail if they are not. Ideal when using any other table to prevent analyzing incorrect columns.''')
##
##    parser.add_argument("-Gp", "--lengthsfilep",    type=str, default=False,
##                        help='''Only used with --bedpe option. Two-column, tab-separated values file with seqnames and lengths for primary seqs. Same as what is called a genome file with bedtools. If BED used rather than faCount output, then this is required. If this given with facount, this supercedes lengths learned from facount.''')
##
##    parser.add_argument("-Sp", "--seqnamecolp",     type=int, default=None,
##                        help='''Only used with --bedpe option. 1-based Column index of --facountp/--seqtablep argument/file to find sequence name. Defaults to --seqnamecol.''')
##    parser.add_argument("-Np", "--ncountcolp",      type=int, default=None,
##                        help='''Only used with --bedpe option. 1-based Column index of --facountp/--seqtablep argument/file to find N-counts in. Defaults to --ncountcol.''')
##    parser.add_argument("-Lp", "--seqlengthcolp",      type=int, default=None,
##                        help='''Only used with --bedpe option. 1-based Column index of --facountp/--seqtablep argument/file to find sequence lengths in. Defaults to --seqlengthcol.''')
##
##    
##    parser.add_argument("-M", "--bedpemapq", type=str, default='0,10',
##                        help='''Only used with --bedpe option.MAPQ cutoffs to analyze BEDPE. Provide comma-separated integers. Default = 0,10''')


    parser.add_argument("-O", "--outputdir", type=str, default=None,
                        help='''Only used with --bedpe option.MAPQ cutoffs to analyze BEDPE. Provide comma-separated integers. Default = 0,10''')

    parser.add_argument("-C", "--cleanup", action='store_true', default=False,
                        help='''Tells you what it is doing. Prints to StdErr.''')
    parser.add_argument("-v", "--verbose", action='store_true', default=False,
                        help='''Tells you what it is doing. Prints to StdErr.''')


    ## Parse Args
    args = parser.parse_args()

    if args.outputdir is None:
        ### RANDOM SEED
        SEED = np.random.randint(10000000,99999999)
        args.outputdir = "assoc_seq_outdir_" + str(SEED)
    


    ## Process Arguments
    ## If providing BED file to specify gaps on query/associated sequences, must also provide lengths file for associated seqs.
    assert not args.bed or (args.bed and args.lengthsfile)
    ## If using seqtable arg, then columns must be specified.
    assert not args.seqtable or (args.seqtable and args.ncountcol is not None and args.seqlengthcol is not None and args.seqnamecol is not None)

    ## if facount table arg used, then columns default to following.
    if args.facount:
        args.ncountcol  = 7
        args.lengthcol  = 2
        args.seqnamecol = 1

    ## if seqtable table arg used, then store that arg in args.facount (see explanation below).
    elif args.seqtable:
        args.facount    = args.seqtable
        ## args.ncountcol and args.lengthcol were required to be explicitly set at command-line.
        ## The code below was written prior to args.seqtable being added.
        ## Thus all that is needed for code to still work is storing the seqtable var in the facount var with diff columns given
        
##    if args.bedpe:
##        ## If providing BED file to specify gaps on target/subject/primary sequences, must also provide lengths file for primary seqs.
##        assert not args.bedp or (args.bedp and args.lengthsfilep)
##        ## Allow the following 3 variables to default to similarly-named variables above if left as None.
##        args.ncountcolp     = args.ncountcol if args.ncountcolp is None else args.ncountcolp
##        args.seqnamecolp    = args.seqnamecol if args.seqnamecolp is None else args.seqnamecolp
##        args.seqlengthcolp  = args.seqlengthcol if args.seqlengthcolp is None else args.seqlengthcolp
##
##        ## If using seqtablep arg, then columns must be specified explicitly or via defaults in previous step.
##        assert not args.seqtablep or (args.seqtablep and args.ncountcolp is not None and args.seqlengthcolp is not None and args.seqnamecolp is not None)
##
##        ## if facountp table arg used, then columns default to following.
##        if args.facountp:
##            args.ncountcolp     = 7
##            args.lengthcolp     = 2
##            args.seqnamecolp    = 1
##
##        ## if seqtablep table arg used, then store that arg in args.facountp (see explanation above).
##        elif args.seqtablep:
##            args.facountp       = args.seqtablep
##            ## See explanation above.

    ## create minscore dictionary
    args.minscore = {}
    args.minscore['pctM'] = args.min_pct_length_in_matches
    args.minscore['pctL'] = args.min_pct_length_in_alignments
    args.minscore['p(M)'] = args.min_prob_given_matches
    args.minscore['p(L)'] = args.min_prob_given_alignments
        
    
    return args




def get_file_lines(fname):
    ''' fname = path to file.
    Note: this def also occurs in partitionOverlappingBEDintervals.py'''
    with open(fname) as fh:
        lines = [e.strip().split() for e in fh.readlines() if not e.startswith('#')]
    return lines

def get_lengths_from_paf(paf, get_target=False, get_query=False):
    lengths = {}
    for aln in paf:
        if get_query:
            lengths[aln[0]] = int(aln[1])
        if get_target:
            lengths[aln[5]] = int(aln[6])
    return lengths


def make_dict_from_table(table, keycol=0, valcol=1, keytype=str, valtype=int):
    '''table = list of uniformly-formatted sublists'''
    tabledict = {}
    for row in table:
        tabledict[ keytype( row[keycol] ) ] = valtype( row[valcol] )
    return tabledict

def make_defaultdict_from_table(table, keycol=0, valcol=1, keytype=str, valtype=int, initialize=False, defaulttype=int):
    '''table = list of uniformly-formatted sublists'''
    tabledict = defaultdict( defaulttype )
    if initialize:
        tabledict = initialize_defaultdict_from_clone(clone         = initialize,
                                                      defaulttype   = defaulttype)
    for row in table:
        tabledict[ keytype( row[keycol] ) ] = valtype( row[valcol] )
    return tabledict


##def make_dict_from_table(table, keycolumn, valcolumn, initialize=False):
##    if initialize:
##        d = initialize_ncount(initialize)
##    for row in table: 
##        d[row[keycolumn]] += int( row[valcolumn] )
##    return d

def get_lengths_from_table(table, namecol=0, lengthcol=1, remove=None):
    '''table = list of uniformly-formatted sublists'''
    tabledict = make_dict_from_table(table   = table,
                                     keycol  = namecol,
                                     valcol  = lengthcol,
                                     keytype = str, 
                                     valtype = int)
    if remove is not None:
        return clean_dict(tabledict, remove)
    return tabledict
            
    


##all_lengths = get_lengths_from_table(table      = get_file_lines( fname = lengths_fname ),
##                                             namecol    = seqnamecol,
##                                             lengthcol  = seqlengthcol)

def clean_dict(d, remove):
    if type(remove) in (str, int, float):
        try:
            d.pop(remove, None)
        except:
            pass
    elif type(remove) in (list, tuple, set):
        for item in remove:
            try:
                d.pop( item , None )
            except:
                pass
    return d


def get_lengths_from_facount(facount, remove_total=True):
    '''facount = list of facount-output-formatted sublists'''
    fadict = get_lengths_from_table(table = facount,
                                    keycol = 0,
                                    valcol = 1,
                                    keytype = str,
                                    valtype = int)
    if remove_total:
        return clean_dict(fadict, 'total') 
    return fadict
    

def initialize_defaultdict_from_clone(clone, defaulttype=int):
    newdict = defaultdict(defaulttype)
    for seqname in clone.keys():
        ## touch key to initialize in defaultdict. ###= 0
        newdict[seqname]        
    return newdict

def initialize_ncount(lengths):
    return initialize_dict_from_clone(clone = lengths)

def get_sum_of_bed_intervals(bed, d):
    for interval in bed:
        d[interval[0]] += int(interval[2]) - int(interval[1])
    return d

def get_ncount_from_bed(bed, lengths):
    ### ...initialize ncount dict
    ncount = initialize_ncount(lengths)
    ### update
    return get_sum_of_bed_intervals(bed, ncount)






def get_ncount_from_table(table, seqnamecol, ncountcol, lengths):
    return make_defaultdict_from_table(table        = table,
                                       keycol       = seqnamecol,
                                       valcol       = ncountcol,
                                       keytype      = str,
                                       valtype      = int, 
                                       initialize   = lengths,
                                       defaulttype  = int)


def update_lengths(lengths, ncount):
    lengths2 = defaultdict(int)
    for seqname in lengths.keys():
        lengths2[seqname] = lengths[seqname] - ncount[seqname]
    return lengths2



def get_gap_adjusted_lengths_from_table(table, namecol, lengthcol, ncountcol, remove=None, initialize=False):
    lengths = make_defaultdict_from_table( table        = table,
                                           keycol       = namecol,
                                           valcol       = lengthcol,
                                           keytype      = str,
                                           valtype      = int,
                                           initialize   = initialize,
                                           defaulttype  = int)
    ncounts = make_defaultdict_from_table( table        = table,
                                           keycol       = namecol,
                                           valcol       = ncountcol,
                                           keytype      = str,
                                           valtype      = int,
                                           initialize   = initialize,
                                           defaulttype  = int)
    if remove is not None:
        lengths = clean_dict(lengths, remove)
        ncounts = clean_dict(ncounts, remove)
    return update_lengths(lengths, ncounts)

def get_gap_adjusted_lengths_from_file(fname, namecol, lengthcol, ncountcol, remove = None):
    return get_gap_adjusted_lengths_from_table(table        = get_file_lines( fname ),
                                               namecol      = namecol,
                                               lengthcol    = lengthcol,
                                               ncountcol    = ncountcol,
                                               remove       = remove) 

def process_paf(paf, bed_fname=False, table_fname=False, get_query=False, get_target=False, seqnamecol=None, seqlengthcol=None, ncountcol=None, lengths_fname=None, verbose=False):
    ''' paf = list object from get_file_lines()
    '''
    verboseMsg(':::    Lengths from PAF.', gate=verbose)

    ### Get lengths of sequences (( 2 ops here ; paf doesnt have all ))
    lengths = get_lengths_from_paf(paf,
                                   get_query = get_query,
                                   get_target = get_target)

    verboseMsg(':::    Gathering all lengths.', gate=verbose)
    if lengths_fname:   ## get lengths from genomefile (uses facount function below b/c both have name:length in 1st and second element positions).
        verboseMsg(':::    - genome file provided.', gate=verbose)
        all_lengths = get_lengths_from_table(table      = get_file_lines( fname = lengths_fname ),
                                             namecol    = 0,
                                             lengthcol  = 1)
    elif table_fname:   ## get lengths from facount table
        verboseMsg(':::    - sequence table file provided.', gate=verbose)
        #all_lengths = get_lengths_from_facount( get_file_lines( fname = table_fname ) )
        all_lengths = get_lengths_from_table(table      = get_file_lines( fname = table_fname ) ,
                                             namecol    = seqnamecol,
                                             lengthcol  = seqlengthcol,
                                             remove     = 'total')


    ### Process BED or FaCount Table
    verboseMsg(':::    Ncounts.', gate=verbose)
    if bed_fname:
        verbose(':::    - BED provided.', gate=verbose)
        nfile = get_file_lines(bed_fname)
        ncount = get_ncount_from_bed(nfile,
                                     lengths = lengths)
    elif table_fname:
        verboseMsg(':::    - Bsequence table file provided.', gate=verbose)
        nfile = get_file_lines(table_fname)
        #ncount = get_ncount_from_table(nfile,
        #                               column = ncol,
        #                               lengths = lengths)
        ncount = get_ncount_from_table(table        = nfile,
                                       seqnamecol   = seqnamecol,
                                       ncountcol    = ncountcol,
                                       lengths      = lengths)
        

    verboseMsg(':::    Update lengths.', gate=verbose)
    ### Get updated lengths of associated sequences
    verboseMsg(':::    - from PAF.', gate=verbose)
    lengths2 = update_lengths(lengths = lengths,
                               ncount = ncount)

    ### Get updated ALL lengths of associated sequences
    verboseMsg(':::    - from ALL.', gate=verbose)
    all_lengths2 = update_lengths(lengths = all_lengths,
                                  ncount = ncount)

    ## RETURN
    return lengths, nfile, ncount, lengths2, all_lengths, all_lengths2


def get_exp_prop_from_lengths(lengths=None, keys=None):
    return normalize_dict_by_sum_vals(d=lengths, keys=keys)

def normalize_dict_by_sum_vals(d=None, keys=None):
    if d is None:
        return None
    subset  = d if keys is None else {key:d[key] for key in keys}
    lsum    = sum(subset.values())
    exp     = {}
    for name in subset.keys():
        if subset[name] == 0:
            exp[name] = subset[name] ## May thwart some 0-Div cases.
        else:
            exp[name] = subset[name]/lsum
    return exp

def weighted_matches_from_partbed_interval(interval):
    '''
    Assumes name, start, end, matches, original_length.

    original_length = length of interval the matches originally came from.
                    = whereas end-start is the length of the partitioned interval
                    = meaning (end-start)/original_length is the proportion of the original interval taken by the partitioned sub-interval.
                    = that allows us to assign a proportional value to matches as well.
    '''
    start       = int(   interval[1] )
    end         = int(   interval[2] )
    matches     = float( interval[3] )
    orig_len    = float( interval[4] )
    return matches*(end - start)/orig_len

def proportion_parental_seq_from_partbed_interval(interval):
    '''
    Assumes name, start, end, parent_sequence_length, __anything__

    __anything__            = Not used. Often is original_length as defined in weighted_matches_from_partbed_interval() is here as dummy var.

    parent_sequence_length  = length of sequence interval is from (i.e. length of named sequence in BED3).

    '''
    start       = int(   interval[1] )
    end         = int(   interval[2] )
    parent_len  = float( interval[3] )
    return (end - start)/parent_len

def get_interval_len_from_bed3(interval):
    '''
    Assumes name, start, end, parent_sequence_length, __anything__

    __anything__            = Not used. Often is original_length as defined in weighted_matches_from_partbed_interval() is here as dummy var.

    parent_sequence_length  = length of sequence interval is from (i.e. length of named sequence in BED3).

    '''
    start       = int(   interval[1] )
    end         = int(   interval[2] )
    return end - start


def create_dict_from_partbed(partbed, function, initialkeys=None):
    '''
    partbed     = list of sublists ; output from run_sliceToBed()
    function    = function to use on sublists in partbed ; should return only one integer or float value that will be added.
    '''
    ## awk 'OFS="\t" {print $1,$2,$3, ($3-$2)/$4}' | \  ### ALT: {print $1,$2,$3,$4*($3-$2)/$5}'
    ##      sortBed -i - | \
    ##      mergeBed -i - -c 4 -o sum | \
    ##      awk 'OFS="\t" {a[$1]+=$4}END{for (e in a) print e, a[e]}' | \
    ##      awk '{gsub(/%/,"\t"); print}' | \
    ##      sort -k1,1 -k3,3nr | \
    ##      tableFilter.py -n 1 -s 3 stdin
    ##
    partbeddict = defaultdict(int)
    if initialkeys is not None:
        for key in initialkeys:
            partbeddict[key] = 0
            
    for interval in partbed:
        name                     = interval[0]
        partbeddict[ name ]     += function(interval)

    return partbeddict


def join_dict_vals_by_list(*args):
    ''' Args should all be dict.
        Can supply any number of dicts.
        Dicts were intended to have str, int, or float values ;
        Also will merge lists or tuples (into lists).
        Dict vals wont work.
        Assumes all dicts share the same keys.
            Asserts as well.
    '''
    d       = defaultdict(list)
    keys    = sorted(list(set([key for arg in args for key in arg.keys()])))

    for arg in args:
        assert set(arg.keys()) == set(keys)
    
    for key in keys:
        for arg in args:
            if type( arg[ key ] ) in (int, float, str, np.float32, np.float64, np.int32, np.int64, np.str_):
                d[key].append( arg[ key ] )
            elif type( arg[ key ] ) == list:
                d[key] += arg[ key ]
            elif type( arg[ key ] ) == tuple:
                d[key] += list(arg[ key ])
            else:
                print("UNEXPECTED TYPE (add to join_dict_vals_by_list):", key, type(key), type(arg[ key ]))
                assert False
            
    return d
    
##def combine_2_partbed_dicts(d1, d2):
##    dout = defaultdict(list)
##    names = list(set(list(d1.keys()) + list(d2.keys())))
##    names.sort()
##    
##    for name in names:
##        dout[name].append( d1[name] )
##        dout[name].append( d2[name] )
##    return dout

def combine_2_partbed_dicts(d1, d2):
    return join_dict_vals_by_list(d1, d2)

def initialize_combined_partbed_dict(initialkeys):
    '''initialkeys = list obj'''
    cpbd = defaultdict(list)
    for key in initialkeys:
        cpbd[key] = []
    return cpbd

def gather_all_seq_pairs_from_paf(paf, minscore=0, scorecol=11):
    ''' paf = paf = list object from get_file_lines() ; list of paf-formatted sublists'''
    l = []
    for aln in paf:
        if float(aln[scorecol]) >= minscore:
            l.append( aln[0] + '%' + aln[5] )
    return l

def initialize_combined_partbed_dict_from_paf(paf, minscore=0, scorecol=11):
    return initialize_combined_partbed_dict(
        gather_all_seq_pairs_from_paf(paf,
                                      minscore,
                                      scorecol)
        )


def get_len_of_lists_in_dict(d):
    return len( d[ d.keys()[0] ] )




def append_combined_partbed_dict(cpbd, pbd):
    ''' for values, can append elements to list and/or merge lists'''
    #A = set(list(cpbd.keys()))
    #B  = set(list(pbd.keys()))
    #A_only = A.difference(B)
    #B_only = B.difference(A)
    #A_and_B = A.intersection(B)
    #list_len = get_len_of_lists_in_dict( cpbd )
    
    ## Address elements that appear in both:
    names = list( cpbd.keys() )
    assert set(names) == set( list(pbd.keys()))
    for name in names:
        if type( pbd[name] ) in (int, float, str):
            cpbd[name].append( pbd[name] )
        elif type( pbd[name] ) == list:
            cpbd[name] += pbd[name]

    ## Address elements that only appear in original dict, cpbd = A only ; these will get 0 values
            
    return cpbd

def normalize_dict(d1, d2, d3=None, pseudo=0):
    '''
    Input:
    d1  = dict with integer or float values to be normalized ; these will be numerator.
    d2  = dict with integer or float values to normalize d1 with; these will be denominator.
        . d2 should have same keys as d1, or d3 is needed for translation.
    d3  = dict with same keys as d1 and values that are keys to d2.
        . default is None
        
    Output:
    d4  = d1 normalized by d2
    '''
    if d3 is None: ## then have d3 return keys for d2
        d3 = {}
        for key in d2.keys():
            d3[key] = key
    d4 = {}
    for key in d1.keys():
        numerator   = float(d1[key])
        denominator = float(d2[ d3[key] ])
        if numerator == 0 and pseudo == 0: ## most division by 0 has 0/0, in which case we can call it 0.
            d4[key] = numerator
        else:
            d4[key] = (float(d1[key]) + pseudo) / (float(d2[ d3[key] ]) + pseudo)
    return d4


def log2ratio_dict(d1, d2, d3=None, pseudo=0):
    return log2_dict_values( normalize_dict(d1, d2, d3, pseudo) )


def median_log2ratio_normalize_dict(d1, d2, d3=None, pseudo=0):
    ## return median-subtracted log2 ratios
    return median_subtract_dict_values( log2ratio_dict(d1, d2, d3, pseudo) )


def log2_dict_values(d):
    return {key:np.log2(d[key]) for key in d.keys()}

def median_subtract_dict_values(d):
    med = np.median( np.array(list(e for e in d.values())) )
    return {key:d[key]-med for key in d.keys()}

def make_translation_dict(d1, delim='%', keepidx=0, inverse=False, listboth=False, otheridx=None):
    '''
    Input:
    d1      = dict with string key names that can be split at a delimiter.
    delim   = delimiter that can split key strings; default = %
    keepidx = index of split key string where piece will be stored as value in new dict
        
    Output:
    d2      = dict with same keys as d1 and values that are keys made by splitting those keys and selecting a piece ot keep.
    '''
    d2 = {}
    for key in d1.keys():
        splitkey    = str(key).split(delim)
        newkey      = splitkey[keepidx]

        if listboth:
            assert otheridx is not None
            assert inverse is False    ## Failing here b/c inverse has no meaning with listboth
            otherkey        = splitkey[otheridx]
            d2[str(key)]    = [newkey, otherkey]
        else:
            if inverse:
                d2[newkey]      = str(key)
            else:
                d2[str(key)]    = newkey
    return d2



def normalize_dict_pipeline(d1, d2, delim='%', keepidx=0):
    '''
    Input:
    d1      = dict with integer or float values to be normalized ; these will be numerator.
    d2      = dict with integer or float values to normalize d1 with; these will be denominator.
                . d2 should have same keys as d1, or d3 is needed for translation.
    delim   = delimiter that can split key strings; default = %
    keepidx = index of split key string where piece will be stored as value in new dict

    Intermediates:
    d3  = dict with same keys as d1 and values that are keys to d2.
    
    Output:
    d4  = d1 normalized by d2
    '''
    d3 = make_translation_dict(d1,
                               delim,
                               keepidx,
                               inverse=False)
    return normalize_dict(d1, d2, d3)



def sum_dict_by_subkeys(d1, delim="%", keepidx=0):
    '''
    Inouts:
    d1  = dict with splittable keys (and target subkeys therein) and integer or float values
    delim   = delimiter that can split key strings; default = %
    keepidx = index of split key string where piece will be stored as value in new dict

    Intermediates:
    d2  = dict with same keys as d1 and values that are subkeys.

    Outputs:
    d3  = dict whose keys are values from d2, which are subkeys from d1 ;
            and values are sums from values in d1 from keys that share the same subkey
    '''
    d2 = make_translation_dict(d1,
                               delim,
                               keepidx,
                               inverse=False)
    d3 = defaultdict(int)
    for key in d1.keys():
        d3[ d2[key] ] += d1[key]
    return d3


def normalize_dict_by_subkey_sums_pipeline(d1, delim='%', keepidx=0):
    '''
    Input:
    d1      = dict with integer or float values to be summed up and normalized
    delim   = delimiter that can split key strings; default = %
    keepidx = index of split key string where piece will be stored as value in new dict

    Intermediates:
    []      = [intermediate] means it is hidden from this function as an intermediate in a function called by this one.
    [d2]    = [sum_dict_by_subkeys] dict with same keys as d1 and values that are subkeys.
    d2      = dict whose keys are values from d2, which are subkeys from d1 ; (so keys are subkeys from d1)
                and values are sums from values in d1 from keys that share the same subkey
    [d3]    = [normalize_dict_pipeline] dict with same keys as d1 and values that are keys to d2.

    Output:
    d4      =  d1 normalized by sums in d2 related by subkeys
    '''
    d2 = sum_dict_by_subkeys(d1         = d1,
                             delim      = delim,
                             keepidx    = keepidx)
    
    return normalize_dict_pipeline(d1       = d1,
                                   d2       = d2,
                                   delim    = delim,
                                   keepidx  = keepidx)


def is_greater(a, b):
    return a > b
def is_lesser(a, b):
    return a < b

def filter_dict(d1, delim='%', keepidx=0, otheridx=1, cmpfxn=is_greater):
    '''
    Input:
    d1          = dict with integer or float values to be compared
    delim       = delimiter that can split key strings; default = %
    keepidx     = index of split key string where piece will be stored as element 0 in list value of new dict
    otheridx    = index of split key string where piece will be stored as element 1 in list value of new dict
    cmpfxn      = [ is_greater | is_lesser ] function that compares new scores to stored scores.

    Intermediate:
    d2

    Output:
    d3
    '''
    d2 = make_translation_dict(d1       = d1,
                               delim    = delim,
                               keepidx  = keepidx,
                               inverse  = False,
                               listboth = True,
                               otheridx = 1)
    d3 = {}

    for key in d1.keys():
        score       = float(d1[key])
        newkey      = d2[key][0]
        otherkey    = d2[key][1]
        try:
            storedscore     = d3[ newkey ][0]
        except:
            d3[ newkey ]    = [ score, otherkey ]
            storedscore     = d3[ newkey ][0]
            
        if cmpfxn(a = score, b = storedscore):
            d3[ newkey ]    = [ score, otherkey ]
            storedscore     = d3[ newkey ][0]
            
    return d3

def assign_seq(d1, allseqs, delim='%', keepidx=0, otheridx=1, cmpfxn=is_greater, minscore=-1):
    '''
    Input:
    d1          = dict with integer or float values to be compared
    delim       = delimiter that can split key strings; default = %
    keepidx     = index of split key string where piece will be stored as element 0 in list value of new dict
    otheridx    = index of split key string where piece will be stored as element 1 in list value of new dict
    cmpfxn      = [ is_greater | is_lesser ] function that compares new scores to stored scores.

    Intermediate:
    d2
    d3
    
    Output:
    d3'          = updates d3 to also report missing seqs as unplaced with score of 1.
    '''
    ## Given scores, select highest
    d3 = filter_dict(d1         = d1,
                     delim      = delim,
                     keepidx    = keepidx,
                     otheridx   = otheridx,
                     cmpfxn     = cmpfxn)
    ## If seq is lower than the minscore, re-assign as unplaced with the complement score (1-score).
    for key in d3.keys():
        if float(d3[ key ][0]) < minscore:
            d3[ key ][0] = 1 - float(d3[ key ][0])
            d3[ key ][1] = 'Unplaced'
    ## For seqs not evaluated, call them Unplaced with probability of 1.
    for key in allseqs:
        try:
            d3[ key ]
        except:
            d3[ key ]    = [ 1.0, 'Unplaced' ]
    return d3
    

def collapse_values(d1, delim=','):
    d2 = {}
    for key in list(d1.keys()):
        d2[key] = (delim).join(str(e) for e in d1[key])
    return d2

def expand_values(d1, delim=','):
    ''' Assumes values stored as delimited string ; converts to list.'''
    d2 = {}
    for key in list(d1.keys()):
        d2[key] = d1[key].split( delim )
    return d2

def subset_values(d1, start=0, end=None):
    ''' for dicts with list or tuple values; returns list or tuple even if len==1.'''
    d2 = {}
    end = end+1 if end == start else end
    for key in list(d1.keys()):
        endval = len(d1[key]) if end is None else end
        d2[key] = d1[key][ start : endval ]
    return d2


def dict_with_lists_to_dict_with_Counter(d1):
    ''' for dicts with list or tuple values; returns list or tuple even if len==1.'''
    d2 = {}
    for key in list(d1.keys()):
        d2[key] = Counter( d1[key] )
    return d2

def normalize_dict_with_Counter(d1):
    ''' takes in dict with Counter() values, and returns dict with elements in counter normalized between 0 and 1.'''
    d2 = {}
    for key in list(d1.keys()):
        denom = float( sum(d1[key].values()) )
        d2[key] = {subkey:d1[key][subkey]/denom for subkey in d1[key].keys()}
    return d2

def dict_to_summary_list(d, maindelim=',', subdelim=':', sortByVals=False):
    if sortByVals:
        #return maindelim.join( [subdelim.join( sorted([str(e) for e in (k,v)], key = lambda x: x[1]) for k,v in d.items())] )
        l = sorted([[str(e) for e in (k,v)] for k,v in d.items()], key = lambda x: x[1], reverse=True)
        return maindelim.join( [subdelim.join(str(e) for e in sl) for sl in l] )
    else:
        return maindelim.join( [subdelim.join(str(e) for e in (k,v)) for k,v in d.items()] )

def summarize_dict_with_Counter(d1):
    ''' takes in dict with Counter() or similar {} values, and returns dict with elements in counter normalized between 0 and 1.'''
    d2 = {}
    for key in list(d1.keys()):
        n = len(d1[key].values())
        d2[key] = [n, dict_to_summary_list(d1[key], maindelim=',', subdelim=':', sortByVals=True)]
    return d2

def expand_values_in_list(d1, delim=','):
    ''' Assumes values stored as list of things including delimited-strings ; converts to list to list of sublists of strings split at delim.'''
    d2 = {}
    for key in list(d1.keys()):
        d2[key] = [ str(e).split(delim) for e in d1[key] ]
    return d2

def subset_values_in_sublists(d1, start, end=None, simplify=False):
    ''' for dicts with list of sublists ; returns dict with lists or subsetted sublists even if a single element. In the latter case, simplify set to True will return just the element.'''
    d2 = {}
    end = end+1 if end == start or end is None else end
    single_element = True if end-start == 1 and simplify else False
    for key in list(d1.keys()):
        if single_element:
            d2[key] = [ l[ start ] for l in d1[key] ]
        else:
            d2[key] = [ l[ start : end ] for l in d1[key] ]
    return d2



def append_dict(d1, add):
    for key in d1.keys():
        d1[key].append(add)
    return d1




def write_out_dict_of_lists(d, outconn=None, delim='\t', collapsevalsdelim=False, valsonly=False):
    if outconn is None:
        outconn = sys.stdout
    names = list(d.keys())
    names.sort()
    for name in names:
        if collapsevalsdelim:
            vals = ','.join(str(e) for e in d[name])
            if valsonly:
                line = (delim).join( str(e) for e in [vals] )
            else:
                line = (delim).join( str(e) for e in [name, vals] )
        else:
            if valsonly:
                line = (delim).join( str(e) for e in d[name])
            else:
                line = (delim).join( str(e) for e in [name] + d[name])
        outconn.write( line + '\n' )
        
def report_combined_partbed_dict(cpbd, outconn, delim='\t', collapsevalsdelim=False):
    return write_out_dict_of_lists(d = cpbd, outconn = outconn, delim = delim, collapsevalsdelim = collapsevalsdelim)



def expanded_assignment_dict_to_reassignment_dict(expanded):
    pass


def get_scores_from_paf_given_mapq(paf, pafmapq, TMPDIR, initialkeys, lengths, minscoredict):
    tmpbed      = TMPDIR + "/MAPQ_" + str(pafmapq) +"-01-tmp-paf.bed"
    tmpslice1   = TMPDIR + "/MAPQ_" + str(pafmapq) +"-02-1-tmp-paf.slice1.txt"
    tmpslice2   = TMPDIR + "/MAPQ_" + str(pafmapq) +"-02-2-tmp-paf.slice2.txt"
    tmppartbed1 = TMPDIR + "/MAPQ_" + str(pafmapq) +"-03-1-tmp-paf.partition1.bed"
    tmppartbed2 = TMPDIR + "/MAPQ_" + str(pafmapq) +"-03-2-tmp-paf.partition2.bed"
    tmpreport1 = TMPDIR + "/MAPQ_" + str(pafmapq) +"-04-tmp-paf.alignments.txt"
    tmpreport2 = TMPDIR + "/MAPQ_" + str(pafmapq) +"-05-tmp-paf.assignments.txt"
    
    ## Make BED-like query-centric file :: as in utilPAF2QueryCentricBEDlike.py -p PAF 
    with open(tmpbed, 'w') as fh:
        bedlikeList = convertPAFtoBEDlike(
            paf         = paf,
            fh          = fh,
            minscore    = pafmapq,
            scorecolumn = 11,
            mode        = 'both')

    ## Convert BED-like file to Slice file, and sort appropriately :: something like :: awk 'OFS="\t" {print $1,$2,"s",NR,$5,$6"\n"$1,$3,"e",NR}'
    slice_list1 = create_sorted_slice_file_from_bed(
        bed_fname   = tmpbed,
        out_fname   = tmpslice1,
        namecol     = 4,
        scorecol    = 5,
        mode        = 'both')
    slice_list2 = create_sorted_slice_file_from_bed(
        bed_fname   = tmpbed,
        out_fname   = tmpslice2,
        namecol     = 18,
        scorecol    = 5,
        mode        = 'both')
    
    ## Create partition BED file (use maxscore method as needed) :: as in partitionOverlappingBEDintervals.py -s SLICEFILE -M
    ## Note: partbed will be None in write mode (on purpose)
    partbed1 = run_sliceToBed(slicefile = tmpslice1,
                   outputfile = tmppartbed1,
                   scoring='max',
                   sortedNames=False,
                   mode='both')
    partbed2 = run_sliceToBed(slicefile = tmpslice2,
                   outputfile = tmppartbed2,
                   scoring='max',
                   sortedNames=False,
                   mode='both')


    sums1 = create_dict_from_partbed( partbed       = partbed1,
                                      function      = weighted_matches_from_partbed_interval,
                                      initialkeys   = initialkeys)
##    sums2 = create_dict_from_partbed( partbed       = partbed2,
##                                      function      = proportion_parental_seq_from_partbed_interval,
##                                      initialkeys   = initialkeys)
    sums2 = create_dict_from_partbed( partbed       = partbed2,
                                      function      = get_interval_len_from_bed3,
                                      initialkeys   = initialkeys)


    sums1norm = normalize_dict_pipeline(d1      = sums1,
                                        d2      = lengths, ##
                                        delim   = '%',
                                        keepidx = 0)
    sums2norm = sums2 ## dummy var for now
    sums2norm = normalize_dict_pipeline(d1      = sums2,
                                        d2      = lengths, ##
                                        delim   = '%',
                                        keepidx = 0)
    
    p1 = normalize_dict_by_subkey_sums_pipeline(d1       = sums1norm,
                                                delim    = '%',
                                                keepidx  = 0)
    p2 = normalize_dict_by_subkey_sums_pipeline(d1       = sums2norm,
                                                delim    = '%',
                                                keepidx  = 0)
    
    pairsums = combine_2_partbed_dicts( d1  = sums1norm,
                                        d2  = sums2norm)
    pairprobs = combine_2_partbed_dicts(d1  = p1,
                                        d2  = p2)

    sub_raw = append_combined_partbed_dict(pairsums, pairprobs)

    with open(tmpreport1,'w') as fh:
        report_combined_partbed_dict(cpbd = sub_raw, outconn = fh)


    ## ASSIGNMENTS
    as1 = assign_seq(d1      = sums1norm,
                    allseqs  = list(lengths.keys()),
                    delim    = '%',
                    keepidx  = 0,
                    otheridx = 1,
                    cmpfxn   = is_greater,
                    minscore = minscoredict['pctM'])
    as2 = assign_seq(d1      = sums2norm,
                    allseqs  = list(lengths.keys()),
                    delim    = '%',
                    keepidx  = 0,
                    otheridx = 1,
                    cmpfxn   = is_greater,
                    minscore = minscoredict['pctL'])
    ap1 = assign_seq(d1      = p1,
                    allseqs  = list(lengths.keys()),
                    delim    = '%',
                    keepidx  = 0,
                    otheridx = 1,
                    cmpfxn   = is_greater,
                    minscore = minscoredict['p(M)'])
    ap2 = assign_seq(d1      = p2,
                    allseqs  = list(lengths.keys()),
                    delim    = '%',
                    keepidx  = 0,
                    otheridx = 1,
                    cmpfxn   = is_greater,
                    minscore = minscoredict['p(L)'])

    as1 = collapse_values( append_dict(as1, add="chrMap,pctM,Q"+str(pafmapq)), delim=",")
    as2 = collapse_values( append_dict(as2, add="chrMap,pctL,Q"+str(pafmapq)), delim=",")
    ap1 = collapse_values( append_dict(ap1, add="chrMap,p(M),Q"+str(pafmapq)), delim=",")
    ap2 = collapse_values( append_dict(ap2, add="chrMap,p(L),Q"+str(pafmapq)), delim=",")

    sub_assign = append_combined_partbed_dict(combine_2_partbed_dicts(as1, as2),
                                     combine_2_partbed_dicts(ap1, ap2))

    with open(tmpreport2,'w') as fh:
        report_combined_partbed_dict(cpbd = sub_assign, outconn = fh)

    return sub_raw, sub_assign



def return_index(x):
    if type(x) in (float, int, str):
        return int(x) - 1
    else:
        ## I expect this to be None, otherwise this is failing silently
        return x 

def verboseMsg(msg, gate=False):
    if gate:
        sys.stderr.write(msg + '\n')





def analyze_paf(paf, args, a_all_lengths2, TMPDIR):
    pafmapqlist = [int(e) for e in args.pafmapq.strip().split(',')]

    
    rawout      = initialize_combined_partbed_dict_from_paf( paf,
                                                        minscore = min(pafmapqlist),
                                                        scorecol = 11)

    assignments = initialize_combined_partbed_dict(initialkeys = list(a_all_lengths2.keys()))
    
    #initialkeys = list(allsums.keys())
    
    for pafmapq in sorted(pafmapqlist):
        sub_raw, sub_assign = get_scores_from_paf_given_mapq( paf          = paf,
                                                              pafmapq      = pafmapq,
                                                              TMPDIR       = TMPDIR,
                                                              initialkeys  = list(rawout.keys()),
                                                              lengths      = a_all_lengths2,
                                                              minscoredict = args.minscore) ## could arguably be a_all_lengths2
        rawout      = append_combined_partbed_dict(rawout, sub_raw)
        assignments = append_combined_partbed_dict(assignments, sub_assign)
        #allsums, pa = append_combined_partbed_dict(allsums, pairsums)
        #report_combined_partbed_dict(cpbd = pa, outconn = sys.stdout)

    ## Pre-TableFilter step in original command-line version. To recreate those final steps: awk '{gsub(/%/,"\t"); print}' | sort -k1,1 -k3,3nr | tableFilter.py -n 1 -s 3 stdin
    #report_combined_partbed_dict(cpbd = allsums, outconn = sys.stdout)

    #report_combined_partbed_dict(cpbd = assignments, outconn = sys.stdout)
    #print(assignments)
##    print(subset_values_in_sublists(
##            d1 = expand_values_in_list(
##                d1 = assignments,
##                delim=','),
##            start = 1,
##            end = 2,
##            simplify=True))
##    quit()

    #print(assignments)
##    print(summarize_dict_with_Counter(
##            normalize_dict_with_Counter(
##                dict_with_lists_to_dict_with_Counter(
##                    d1 = subset_values_in_sublists(
##                        d1 = expand_values_in_list(
##                            d1 = assignments,
##                            delim=','),
##                        start = 1,
##                        end = 2,
##                        simplify=True),
##                        )
##                    )
##                )
##          )


    expanded = expand_values_in_list(
                            d1 = assignments,
                            delim=',')
    
    subsetted = subset_values_in_sublists(
                        d1 = expanded,
                        start = 1,
                        end = 2,
                        simplify=True)

    
    normalizedCounterDict = normalize_dict_with_Counter( dict_with_lists_to_dict_with_Counter( d1 = subsetted ) )
    summarizedCounterDict = summarize_dict_with_Counter( normalizedCounterDict )



    finalreport = TMPDIR + "/Final-report-paf.assignment-summary.txt"
    with open(finalreport, 'w') as fh:
        report_combined_partbed_dict(
            summarizedCounterDict,
            outconn = fh)

    #"--min_mode_of_assignment"
    #print(summarizedCounterDict)
    final_assignments = [[key,summarizedCounterDict[key][1].split(':')[0]] if float(summarizedCounterDict[key][1].split(':')[1].split(',')[0]) > args.min_mode_of_assignment else [key, 'Unplaced'] for key in summarizedCounterDict.keys() ]
    return final_assignments
    

        


    
def is_gzipped(fname):
    return fname.endswith('.gz')

def get_file_connection(fname, mode):
    if is_gzipped( fname ):
        return io.TextIOWrapper(gzip.GzipFile( fname, mode ))
    else:
        return open( fname, mode )


def merge_dicts(x, y=None):
    if y is None:
        return x
    ## Else try merges:
    try:
        return x | y            ## python >= 3.9
    except:
        pass
    try:
        return {**a, **b}      ## python >= 3.4
    except:
        pass
    try:
        x.update(y)      ## catchall
        return x
    except:
        pass
    return None
    


def main():
    args = parse_args()

    os.mkdir( args.outputdir )
    
    ### Get PAF
    verboseMsg(':::    Reading Paf.', gate=args.verbose)
    paf = get_file_lines(fname = args.paf)
    

    ## PROCESS
    verboseMsg(':::    Processing Paf.', gate=args.verbose)
    alengths, nfile, ncount, alengths2, a_all_lengths, a_all_lengths2 = process_paf(paf = paf,
                                                                                    bed_fname       = args.bed,
                                                                                    table_fname     = args.facount,
                                                                                    get_query       = True,
                                                                                    get_target      = False,
                                                                                    seqnamecol      = return_index( args.seqnamecol ),
                                                                                    seqlengthcol    = return_index( args.seqlengthcol ),
                                                                                    ncountcol       = return_index( args.ncountcol ),
                                                                                    lengths_fname   = args.lengthsfile,
                                                                                    verbose         = args.verbose)
    verboseMsg(':::    Done processing Paf.', gate=args.verbose)


    ### Get BEDPE
##    if args.bedpe:
##        verboseMsg(':::    Reading Bedpe.', gate=args.verbose)
##        ### Get bedpe lines
##        bedpe = get_file_lines(args.bedpe)
##
##        ### Get information to help compute expected proportions of links
##        plengths, nfilep, ncountp, plengths2, p_all_lengths, p_all_lengths2 = process_paf(paf_fname     = paf,
##                                                                                          bed_fname     = args.bedp,
##                                                                                          table_fname   = args.facountp,
##                                                                                          get_query     = False,
##                                                                                          get_target    = True,
##                                                                                          seqnamecol    = return_index( args.seqnamecolp ),
##                                                                                          seqlengthcol  = return_index( args.seqlengthcolp ),
##                                                                                          ncountcol     = return_index( args.ncountcolp ),
##                                                                                          lengths_fname = args.lengthsfilep)
##
##
##
##            
##        ### Get expected proportions
##        exp_link_props = get_exp_prop_from_lengths( p_all_lengths2 )
        


    ### RANDOM SEED
    #SEED = np.random.randint(10000000,99999999)
    #TMPDIR = "assoc_seq_tmpdir_" + str(SEED)
    #os.mkdir(TMPDIR)


    ##############################################################################
    ''' ANALYZE PAF '''
    ##############################################################################
    ### CREATE BED FROM PAF
        ## awk '$12 > Q {OFS="\t"; print $1"%"$6, $3,$4, $2}' PAF | sortBed -i - | mergeBed -i - -c 4 -o distinct | awk 'OFS="\t" {a[$1"%"$4]+=($3-$2)}END{ for (e in a) print e,a[e]}' | awk 'OFS="\t" {gsub(/%/,"\t"); print}' | sort -k1,1 -k4,4nr
        ## utilPAF2QueryCentricBEDlike.py -p del1.paf | awk 'OFS="\t" {print $1,$2,"s",NR,$5,$6"\n"$1,$3,"e",NR}' | sort -k1,1 -k2,2n | partitionOverlappingBEDintervals.py -s - -M | awk 'OFS="\t" {print $1,$2,$3,$4*($3-$2)/$5}' | sortBed -i - | mergeBed -i - -c 4 -o sum | awk 'OFS="\t" {a[$1]+=$4}END{for (e in a) print e, a[e]}' | awk '{gsub(/%/,"\t"); print}' | sort -k1,1 -k3,3nr | tableFilter.py -n 1 -s 3 stdin > res1
        ## utilPAF2QueryCentricBEDlike.py -p del1.paf | awk 'OFS="\t" {print $1,$2,"s",NR,$19,$6"\n"$1,$3,"e",NR}' | sort -k1,1 -k2,2n | partitionOverlappingBEDintervals.py -s - -M | awk 'OFS="\t" {print $1,$2,$3, ($3-$2)/$4}' | sortBed -i - | mergeBed -i - -c 4 -o sum | awk 'OFS="\t" {a[$1]+=$4}END{for (e in a) print e, a[e]}' | awk '{gsub(/%/,"\t"); print}' | sort -k1,1 -k3,3nr | tableFilter.py -n 1 -s 3 stdin > res2

    final_assignments = analyze_paf(paf, args, a_all_lengths2, args.outputdir)

    for f in final_assignments:
        f.append( a_all_lengths[f[0]] )

    final_assignments = sorted(final_assignments, key = lambda x: (x[1], -int(x[2])) )
    
    with open(args.outputdir + "/Final-assignments.txt", 'w') as fh:
        for f in final_assignments:
            fh.write( '\t'.join(str(e) for e in f) + '\n' )

    outtypes = sorted(list(set([f[1] for f in final_assignments])))
    d = {}
    for e in outtypes:
        d[e] = open(args.outputdir + "/Final-assignments." + e +'.txt', 'w')
    d['Placed'] = open(args.outputdir + "/Final-assignments." + 'Placed' +'.txt', 'w')

    for f in final_assignments:
        d[ f[1] ].write( f[0] + '\n' )
        if f[1] != 'Unplaced':
            d[ 'Placed' ].write( f[0] + '\n' )

    for key in d.keys():
        d[key].close()



    ##############################################################################
    ''' ANALYZE BEDPE '''
    ##############################################################################



    ##############################################################################
    ''' SYNTHESIZE OUTPUT '''
    ##############################################################################



    ##############################################################################
    ''' CLEAN UP '''
    ##############################################################################
    #if args.cleanup:
    #    shutil.rmtree(TMPDIR)


##############################################################################
##############################################################################
##############################################################################
'''
EXECUTE
'''
##############################################################################
##############################################################################
##############################################################################
if __name__ == "__main__":
    main()
