#!/usr/bin/env python3

import argparse, sys, os, shutil, pybedtools, gzip
from collections import defaultdict, Counter
import numpy as np

## SCIARA TOOLS IMPORTS
SCITOOLS_BIN_PATH = os.path.dirname(__file__)
if not SCITOOLS_BIN_PATH in sys.path:
    sys.path.append(SCITOOLS_BIN_PATH)
from utilPAF2QueryCentricBEDlike import convertPAFtoBEDlike
from partitionOverlappingBEDintervals import SliceLine, slice3, run_sliceToBed, sortbed, create_sorted_slice_file_from_bed
from associateSequences import *

##############################################################################
''' FUNCTIONS '''
##############################################################################

def parse_args():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter, description='''
    pass        
    ''')

    parser.add_argument("-b", "--bedpe", type=str, required=True, help='''Required bedpe file (bamToBed output style).''')
    
    parser.add_argument("-a", "--associaseqtable", "--facount",   type=str, default=False,
                         help='''Tab-sep text file with seq names, lengths, and N-counts for associated seqs. Assumes faCount output, seq name in column 1, seq lengths in column 2, and N count in column 7. Use --seqnamecol, --seqlengthcol, and --ncountcol to change.''')
    parser.add_argument("-p", "--primaryseqtable", "--facountp",   required=True, type=str, default=False,
                         help='''Tab-sep text file with seq names, lengths, and N-counts for associated seqs. Uses same columns used for --assocseqtable.''')

    parser.add_argument("-S", "--seqnamecol",   type=int, default=1,
                        help='''1-based Column index of --facount/--seqtable argument/file to find sequence name.''')
    parser.add_argument("-L", "--seqlengthcol",    type=int, default=2,
                        help='''1-based Column index of --facount/--seqtable argument/file to find sequence lengths in.''')
    parser.add_argument("-N", "--ncountcol",    type=int, default=7,
                        help='''1-based Column index of --facount/--seqtable argument/file to find N-counts in.''')



    parser.add_argument("-m", "--minmapq", type=float, default=float('-inf'),
                        help='''Minimum MAPQ cutoff. Default = -inf.''')
    parser.add_argument("-M", "--maxmapq", type=float, default=float('inf'),
                        help='''Maximum MAPQ cutoff. Default = inf.''')
    parser.add_argument("--pseudo", type=float, default=1,
                        help='''Pseudocounts for links between 2 sequences. Float. Default = 1.0''')

    parser.add_argument("-O", "--outdir", "--outputdir", type=str, default=None,
                        help='''Only used with --bedpe option.MAPQ cutoffs to analyze BEDPE. Provide comma-separated integers. Default = 0,10''')


    parser.add_argument("-C", "--cleanup", action='store_true', default=False,
                        help='''Not meaningful at the moment.''')
    parser.add_argument("-v", "--verbose", action='store_true', default=False,
                        help='''Tells you what it is doing. Prints to StdErr.''')


    ## Parse Args
    args = parser.parse_args()


    if args.outdir is None:
        ### RANDOM SEED
        SEED = np.random.randint(10000000,99999999)
        args.outputdir = "assoc_seq_by_links_outdir_" + str(SEED)
        
    return args


class HiC_BEDPE_Dict(object):
    def __init__(self, fname):
        '''
        fname       = path to BEDPE of HIC links 
                        a = {assocseq:{pseq:{mapq0:10}}}
        
        '''
        self.fname      = fname
        self.hic        = {}        ## keys are seq pair lists [X, Y], values are
        self.pairs      = defaultdict(int)
        self.totalpairs = 0
        self.seq        = {}
        self.reads      = defaultdict(int)
        self.totalreads = 0
        self._hic_bedpe_to_dict()
        

    def _hic_bedpe_to_dict(self):
        hic = get_file_connection(fname = self.fname, mode = 'r')
        for line in hic:
            ## Process line
            line        = [str(e) for e in line.strip().split()] ## str(e) helps with gzip files.
            seqpair     = sorted( [line[0], line[3]] )
            mapq        = int(line[7])

            ## Add
            self._add_to_hic_dict( seqpair    , mapq , add = 1 )
            self._add_to_seq_dict( seqpair[0] , mapq , add = 1 )
            self._add_to_seq_dict( seqpair[1] , mapq , add = 1 )
        hic.close()

    def _add_to_hic_dict(self, seqpair, mapq, add=1):
        seqpair = sorted(seqpair)
        mapq    = int(mapq)
        add     = int(add)
        ## Add to HiC pair counts
        try:
            self.hic[seqpair[0]]
        except:
            self.hic[seqpair[0]] = {}
        try:
            self.hic[seqpair[0]][seqpair[1]]
        except:
            self.hic[seqpair[0]][seqpair[1]]    =   defaultdict(int)
        #try:
        #    self.hic[seqpair[0]][seqpair[1]][mapq] += 1
        #except:
        #    self.hic[seqpair[0]][seqpair[1]][mapq] = 1
        self.hic[seqpair[0]][seqpair[1]][mapq]  +=  add

        ## Add to total number pairs at given mapq level 
        self.pairs[mapq]    += add
        
        
        ## Add to Total count
        self.totalpairs     +=  add


    def _add_to_seq_dict(self, seq, mapq, add=1):
        mapq    = int(mapq)
        add     = int(add)
        ## Add to individual seq counts
        try:
            self.seq[seq]
        except:
            self.seq[seq] = defaultdict(int)
            
        self.seq[seq][mapq] +=  add

        ## Add to total number pairs at given mapq level 
        self.reads[mapq]    += add

        
        ## Add to Total count
        self.totalreads     +=  add


    ########################################################################################################################
    ''' SE / READ ANALYSIS METHODS '''
    ########################################################################################################################
    def get_total_reads(self, minmapq = float('-inf'), maxmapq = float('inf')):
        return sum([self.reads[mapq] for mapq in list( self.reads.keys() ) if mapq >= minmapq and mapq <= maxmapq] )
    
    def get_seq_count(self, seq, mapq):
        ## Unencountered seq,mapq pairs automatically return 0
        return self.seq[seq][mapq]

    def _try_to_get_all_mapq_values_seen_for_seq(self, seq):
        try:
            return self.seq[seq].keys()
        except:
            return []
        
    def get_mapq_values_seen_for_seq(self, seq, minmapq = float('-inf'), maxmapq = float('inf')):
        return sorted([e for e in list(self._try_to_get_all_mapq_values_seen_for_seq(seq)) if e >= minmapq and e <= maxmapq] )

    def get_seq_count_list_given_mapq_thresholds(self, seq, minmapq = float('-inf'), maxmapq = float('inf')):
        return [ self.get_seq_count(seq, mapq) for mapq in self.get_mapq_values_seen_for_seq(seq,
                                                                                             minmapq,
                                                                                             maxmapq) ]
    def get_sum_seq_counts_given_mapq_thresholds(self, seq, minmapq = float('-inf'), maxmapq = float('inf')):
        return sum(
            self.get_seq_count_list_given_mapq_thresholds(
                seq,
                minmapq,
                maxmapq)
            )
    def get_seq_count_dict_given_mapq_thresholds(self, seqs, minmapq = float('-inf'), maxmapq = float('inf')):
        '''
        Input:
        seqs        = list of sequence names to get summed counts for given min and max mapq values
        minmapq     = minimum mapq
        maxmapq     = maximum mapq
        '''
        return {seq:self.get_sum_seq_counts_given_mapq_thresholds(seq,
                                                                  minmapq,
                                                                  maxmapq) for seq in seqs}
    def get_seq_count_proportion_dict_given_mapq_thresholds(self, seqs, minmapq = float('-inf'), maxmapq = float('inf')):
        '''
        Input:
        seqs        = list of sequence names to get summed counts for given min and max mapq values
        minmapq     = minimum mapq
        maxmapq     = maximum mapq

        Uses normalize_dict_by_sum_vals(), a function written outside of this class.
        '''
        return  normalize_dict_by_sum_vals(
            d       = self.get_seq_count_dict_given_mapq_thresholds(
                        seqs,
                        minmapq,
                        maxmapq),
            keys    = None) ## Keys = None here just means to not further subset this list.


        
    ########################################################################################################################
    ''' PE / LINK / PAIR ANALYSIS METHODS '''
    ########################################################################################################################
    def get_total_links(self, minmapq = float('-inf'), maxmapq = float('inf')):
        return sum([self.pairs[mapq] for mapq in list( self.pairs.keys() ) if mapq >= minmapq and mapq <= maxmapq] )
    
    def get_link_count(self, seq1, seq2, mapq):
        seqpair = sorted([seq1, seq2])
        try:
            return self.hic[seqpair[0]][seqpair[1]][mapq]
        except:
            ## This pair and mapq was not encountered ; add it with 0 count, and return.
            self._add_to_hic_dict(self, seqpair, mapq, add=0)
            return self.hic[seqpair[0]][seqpair[1]][mapq]

    def _try_to_get_all_mapq_values_seen_for_seqpair(self, seqpair):
        try:
            return self.hic[seqpair[0]][seqpair[1]].keys()
        except:
            return []

    def get_mapq_values_seen_for_seqpair(self, seq1, seq2, minmapq = float('-inf'), maxmapq = float('inf')):
        seqpair = sorted([seq1, seq2])
        return sorted([e for e in list( self._try_to_get_all_mapq_values_seen_for_seqpair(seqpair) ) if e >= minmapq and e <= maxmapq] )

    def get_link_count_list_for_seqpair_given_mapq_thresholds(self, seq1, seq2, minmapq = float('-inf'), maxmapq = float('inf')):
        seqpair = sorted([seq1, seq2]) ## This is redundant, but just for security.
        return [ self.get_link_count(seq1, seq2, mapq) for mapq in self.get_mapq_values_seen_for_seqpair(seqpair[0],
                                                                                                         seqpair[1],
                                                                                                         minmapq,
                                                                                                         maxmapq) ]
    def get_sum_link_counts_for_seqpair_given_mapq_thresholds(self, seq1, seq2, minmapq = float('-inf'), maxmapq = float('inf')):
        seqpair = sorted([seq1, seq2])
        return sum(
            self.get_link_count_list_for_seqpair_given_mapq_thresholds(
                seqpair[0],
                seqpair[1],
                minmapq,
                maxmapq)
            )
    def get_link_count_dict_for_seq_given_mapq_thresholds_and_seqlist(self, seq1, seqs, minmapq = float('-inf'), maxmapq = float('inf')):
        '''
        Input:
        seq1        = reference seq to get link counts for...
        seqs        = list of sequence names to get summed counts for given min and max mapq values
        minmapq     = minimum mapq
        maxmapq     = maximum mapq
        '''
        return {seq2:self.get_sum_link_counts_for_seqpair_given_mapq_thresholds(seq1,
                                                                                seq2,
                                                                                minmapq,
                                                                                maxmapq) for seq2 in seqs}
    def get_link_count_proportion_dict_for_seq_given_mapq_thresholds_and_seqlist(self, seq1, seqs, minmapq = float('-inf'), maxmapq = float('inf')):
        '''
        Input:
        seq1         = reference seq to get link counts for...
        seqs        = list of sequence names to get summed counts for given min and max mapq values
        minmapq     = minimum mapq
        maxmapq     = maximum mapq

        Uses normalize_dict_by_sum_vals(), a function written outside of this class.
        '''
        return  normalize_dict_by_sum_vals(
            d       = self.get_link_count_dict_for_seq_given_mapq_thresholds_and_seqlist(
                        seq1,
                        seqs,
                        minmapq,
                        maxmapq),
            keys    = None)



class HiC_Analysis_Given_MAPQ_thresholds(object):
    def __init__(self, bedpedict, pseqlength, aseqlength=None, minmapq = float('-inf'), maxmapq = float('inf'), pseudo=1, outdir=None):
        ''' bedpedict       =   HiC_BEDPE_Dict object
            [a|p]seqlength  =   dict with seqname keys and length values ;
                                - p = primary; associated = associated 
                                    - P and A should not share any keys
                                    - giving both P and A will yield further analysis associating each seq in A with 1 seq from P
                                    - only use P with all seqs if not interested in that analysis.
                                - assumes all seqs represented in bedpedict are also here
                                - recommend using modified seqlenghts = actual_seq_len - n_counts
                                    - this modification need be done prior to HiC_Analysis()

        Depends directly on:
            - merge_dicts
            - get_exp_prop_from_lengths
            - HiC_BEDPE_Dict
            - median_ratio_normalize_dict
            - and indirectly on any dependencies for those functions and classes
        '''
        self.bedpedict              = bedpedict
        self.pseqs                  = sorted( list( pseqlength.keys() ) )
        self.aseqs                  = sorted( list( aseqlength.keys() ) ) if aseqlength is not None else None
        self.allseqs                = sorted( self.pseqs + self.aseqs )
        self.minmapq                = minmapq
        self.maxmapq                = maxmapq
        self.pseudo                 = pseudo
        self.outdir                 = outdir
        self.totallinks             = self.bedpedict.get_total_links(
                                            minmapq  = self.minmapq,
                                            maxmapq  = self.maxmapq )
        self.totalreads             = self.bedpedict.get_total_reads(
                                            minmapq  = self.minmapq,
                                            maxmapq  = self.maxmapq )
        if self.aseqs is not None:
            ## Ensure aseqs and pseqs are non-overlapping sets.
            pseqs_and_aseqs_are_nonoverlapping  = len(set( self.pseqs ).intersection(set( self.aseqs ))) == 0
            #assert pseqs_and_aseqs_are_nonoverlapping

        self.seqlength                  = merge_dicts( x = pseqlength, y = aseqlength )
        
        self.expected_proportions       = get_exp_prop_from_lengths(lengths = self.seqlength,
                                                                    keys    = None)

        self.expected_counts            = self._get_expected_seq_count_proportion_dict_given_mapq_thresholds()

        self.expected_proportions_P     = get_exp_prop_from_lengths(lengths = pseqlength,
                                                                    keys    = None)

        self.expected_proportions_A     = get_exp_prop_from_lengths(lengths = aseqlength,
                                                                    keys    = None)

        self.observed_counts            = self.bedpedict.get_seq_count_dict_given_mapq_thresholds(
                                                seqs    = list( self.expected_proportions.keys() ),
                                                minmapq = self.minmapq,
                                                maxmapq = self.maxmapq)
        
        self.observed_proportions       = self.bedpedict.get_seq_count_proportion_dict_given_mapq_thresholds(
                                                seqs    = list( self.expected_proportions.keys() ),
                                                minmapq = self.minmapq,
                                                maxmapq = self.maxmapq)
        

        self.observed_proportions_P     = None

        self.observed_proportions_A     = None

        self.obs_to_exp_ratio                       = normalize_dict(
                                                            d1      = self.observed_counts,
                                                            d2      = self.expected_counts,
                                                            d3      = None,
                                                            pseudo  = self.pseudo)

        self.log2_obs_to_exp_ratio                  = log2ratio_dict(
                                                            d1      = self.observed_counts,
                                                            d2      = self.expected_counts,
                                                            d3      = None,
                                                            pseudo  = self.pseudo)
        self.mednorm_log2_obs_to_exp_ratio           = median_log2ratio_normalize_dict(
                                                            d1      = self.observed_counts,
                                                            d2      = self.expected_counts,
                                                            d3      = None,
                                                            pseudo  = self.pseudo)



        if self.outdir is not None:
            self.covfile    = self.outdir + '/' + '01-coverage-analysis.txt'
            self.linkfile   = self.outdir + '/' + '02-linkage-analysis.txt'
        
            with open(self.covfile, 'w') as fh:
                self.report_cov_analysis( outconn = fh )

            with open(self.linkfile, 'w') as fh:
                self.report_link_analysis( outconn = fh )
        else:
            self.report_cov_analysis( outconn = sys.stdout )
            self.report_link_analysis( outconn = sys.stdout )
            
        
        ##LINKS -- look at all
        ##self.expected_link_proportions_with_intra   = 00000####

        ## Get observed proportions of reads for each sequence
        ## Analyze obs vs exp w/r/t all seqs and the A and/or P subsets
        ## Inter-chromo probs

        ## Expected number of links between A and B is the independent joint prob: p(A)*p(B) = p(B|A)*p(A) = p(A|B)*p(B)
        ##      not just marginalized p(B|A) or p(A|B) ; and not just p(A) or p(B)
        ##      also p(A) and p(B) is computed by independent coverage of each ; not by length of each
        ##      So the expected joint proportion for each is: obscovprop(A)*obscovprop(B)
        ##      The observed joint proportion = #Links(A,B)/#Links(All)
        ##      Can also get exp and obs marginalized to just A and B by:
        ##      Exp(A,B) = p(A^B)/p(AUB) = p(A)*p(B)/(p(A)+p(B)-p(A)*p(B))
        ##      Obs(A,B) = #Links(A,B)/#Links(A or B).. latter is like jaccard..
        ##      The two approaches are related by the following:
        ##          p(A) = #Links_w_A / #Links
        ##          p(B) = #Links_w_B / #Links
        ##          E(A^B) = p(A)*p(B) = (#Links_w_A * #Links_w_B) / (#Links * #Links)
        ##          E(AUB) = p(A) + p(B) - p(A)*p(B)
        ##                 = (#Links/#Links)*[(#Links_w_A + #Links_w_B)/#Links] - [(#Links_w_A * #Links_w_B) / (#Links * #Links)]
        ##                 = [(#Links*#Links_w_A + #Links*#Links_w_B)/(#Links * #Links)] - ""
        ##                 = (#Links*#Links_w_A + #Links*#Links_w_B - #Links_w_A*#Links_w_B)/ (#Links * #Links)
        ##          O(A^B) = #Links_w_A_and_B / #Links
        ##          
        ##          O(AUB) = (#Links_w_A + #Links_w_B - #Links_w_A_and_B) / #Links
        ##
        ##          Exp(Jaccard)    = E(A^B)/E(AUB)
        ##                          = [(#Links_w_A * #Links_w_B) / (#Links * #Links)] / [(#Links*#Links_w_A + #Links*#Links_w_B - #Links_w_A*#Links_w_B) / (#Links * #Links)]
        ##                          = [(#Links_w_A * #Links_w_B) / (#Links * #Links)] * [(#Links * #Links) / (#Links*#Links_w_A + #Links*#Links_w_B - #Links_w_A*#Links_w_B)]
        ##                                  cancels (#Links * #Links)/(#Links * #Links)
        ##                          = (#Links_w_A * #Links_w_B) / (#Links*#Links_w_A + #Links*#Links_w_B - #Links_w_A*#Links_w_B)
        ##          Obs(Jaccard)    = [#Links_w_A_and_B / #Links] / [(#Links_w_A + #Links_w_B - #Links_w_A_and_B) / #Links]
        ##                          = [#Links_w_A_and_B / #Links] * [#Links / (#Links_w_A + #Links_w_B - #Links_w_A_and_B)]
        ##                                  cancels #Links/#Links
        ##                          = #Links_w_A_and_B /  (#Links_w_A + #Links_w_B - #Links_w_A_and_B)

        ##  above is mostly correct...
        ##  some issues arose:
        ##      E(A^B) is from read proportions ; there are 2x more reads than pairs
        ##      so either:
        ##          O(A^B) needs to be: #Links_w_A_and_B / #Reads / #Reads = ......
        ##              I guess the way to think of it is #Reads_w_A_and_B / #Reads ... or #Links_w_A_and_B / 2*#Links
        ##      Or:
        ##          E(A^B) needs to be 2*p(A)*p(B)

        ## ACTUALLY
        ## We don't need the joint prob A^B... we need the conditional prob of observing B given A.
        ##  We know we have A, so we don't need p(A)
        ##  If it is all random and indpendent, then p(B|A) = p(B)
        ## So the expected proportions are just the observed cov proportions
        ## The observed p(B|A) = #A^B / #A
        ## and Im not trying to do O(joint)/E(joint) --- Im trying to do O(B|A)/E(B|A)
    
        ##LINKS -- look at all
        ##self.expected_link_proportions_with_intra   = 00000####

        ## Due to lack of dev time, for now I am only doing the primary vs assoc link analysis here.

        ## For each assoc_seq, get exp and obs proportions marginalized to only the primary seqs.
        ##      In the future, you can make an all-x-all interaction rate table .... then one can do clustering analyses as well
        ##      ...thus some assoc tigs not well-assigned might clister with those that are.

        
        ##
    
        

            
    def report_cov_analysis(self, outconn=None, delim='\t', collapsevalsdelim=False):
        write_out_dict_of_lists(d                   = join_dict_vals_by_list(self._get_total_read_count_dict(),
                                                                             self.expected_counts,
                                                                             self.observed_counts,
                                                                             self.expected_proportions,
                                                                             self.observed_proportions,
                                                                             self.obs_to_exp_ratio,
                                                                             self.log2_obs_to_exp_ratio,
                                                                             self.mednorm_log2_obs_to_exp_ratio),
                                outconn             = outconn,
                                delim               = delim,
                                collapsevalsdelim   = collapsevalsdelim)


    def report_link_analysis(self, outconn=None):
        de = {}
        do = {}

        
        for seq1 in sorted(self.allseqs):
            do[seq1] = self.bedpedict.get_link_count_dict_for_seq_given_mapq_thresholds_and_seqlist(seq1,
                                                                                      seqs    = self.allseqs,
                                                                                      minmapq = self.minmapq,
                                                                                      maxmapq = self.maxmapq)


        
        for seq1 in sorted(self.aseqs):
            d = {}
            o = {}
            e = {}
            for seq2 in sorted(self.pseqs):
                exp_prop_seq2  = self.observed_proportions[seq2]
                obs_nseq1 = sum(do[seq1].values())
                exp_nseq2 = round( exp_prop_seq2*obs_nseq1, 3 )
                obs_joint = do[seq1][seq2]
                obs_nseq2 = obs_joint ## dummy/redundant, but like the naming for some things
                obs_prop_seq2  = ((obs_joint+self.pseudo) / (obs_nseq1+self.pseudo)) if obs_nseq1 == 0 else obs_joint / obs_nseq1
                #ratio1 = obs_joint / (exp_prop_seq2*obs_nseq1)
                #ratio = obs_prop_seq2 / exp_prop_seq2
                d[seq2] = [seq1, seq2, obs_nseq1, exp_nseq2, obs_joint, exp_prop_seq2, obs_prop_seq2]
                o[seq2] = obs_nseq2
                e[seq2] = exp_nseq2

            ratio = normalize_dict(d1    = o,
                                   d2    = e,
                                   d3    = None,
                                   pseudo = self.pseudo)


            write_out_dict_of_lists(d = join_dict_vals_by_list(d,
                                       ratio,
                                       log2_dict_values(ratio),
                                       median_subtract_dict_values(log2_dict_values(ratio))),
                                    outconn = outconn,
                                    delim = '\t',
                                    collapsevalsdelim   = False,
                                    valsonly            = True)
            
        #print( '\t'.join(str(e) for e in l) )
        #print(do)

    def _get_expected_seq_count_proportion_dict_given_mapq_thresholds(self):
        '''
        '''
        return {key:round(self.expected_proportions[key] * self.totalreads, 2) for key in self.expected_proportions.keys()}

    def _get_total_read_count_dict(self):
        '''
        '''
        return {key:self.totalreads for key in self.expected_proportions.keys()}
        
        

def main():
    args = parse_args()
    if args.outdir is not None:
        os.mkdir( args.outdir )

    ### 
    verboseMsg(':::    Constructing links object.', gate=args.verbose)
    bedpedict = HiC_BEDPE_Dict( fname = args.bedpe )

    

    ## PROCESS
    verboseMsg(':::    Preparing other inputs for link analysis.', gate=args.verbose)

    ### Get information to help compute expected proportions of links
    primary_lengths = get_gap_adjusted_lengths_from_file(fname      = args.primaryseqtable,
                                                         namecol    = args.seqnamecol - 1,
                                                         lengthcol  = args.seqlengthcol - 1,
                                                         ncountcol  = args.ncountcol - 1,
                                                         remove     = 'total' )
    associa_lengths = get_gap_adjusted_lengths_from_file(fname      = args.associaseqtable,
                                                         namecol    = args.seqnamecol - 1,
                                                         lengthcol  = args.seqlengthcol - 1,
                                                         ncountcol  = args.ncountcol - 1,
                                                         remove     = 'total' )
    
    verboseMsg(':::    Analyzing links.', gate=args.verbose)
    link_analysis   = HiC_Analysis_Given_MAPQ_thresholds( bedpedict     = bedpedict,
                                                          pseqlength    = primary_lengths,
                                                          aseqlength    = associa_lengths,
                                                          minmapq       = args.minmapq,
                                                          maxmapq       = args.maxmapq,
                                                          pseudo        = args.pseudo,
                                                          outdir        = args.outdir)
    
    verboseMsg(':::    Done.', gate=args.verbose)
    
    


        
    ##############################################################################
    ''' CLEAN UP '''
    ##############################################################################
    if args.cleanup:
        shutil.rmtree(TMPDIR)













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

