#!/usr/bin/env python3

import argparse, sys

def parse_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter, description='''
        INPUT
            PAF

        OUTPUT
            BEDLIKE (BED3 coords followed by various annotation columns)
                1   name            = sequence name
                2   start           = 0-based start coord
                3   end             = end coord
                4   mapq            = mapq
                5   matches         = number matches reported in column 10 of PAF
                6   length          = end - start
                7   alnLen          = length of alignment reported in column 11 of PAF
                8   pident          = matches / alnLen
                9   pdiver          = 1 - pident
                10  pmatch          = matches / length
                11  ambig           = number ambiguous bases in the alignment (nn tag in PAF with mimimap2 -c )
                12  lengthUpdated   = length - ambig
                13  alnLenUpdated   = alnLen - ambig
                14  pident2         = matches / alnLenUpdated
                15  pdiver2         = 1 - pident2
                16  pmatch2         = matches / lengthUpdated
                17  dv              = divergence reported by Minimap2; de if present, dv otherwise.
                18  cv              = 1 - dv
                19  parent_seqlen   = length of sequence interval comes from (eg. col 2 in PAF)
        
    ''')
    parser.add_argument("-p", "--paf", type=str, required=True, help='''PAF alignment file of assoc seqs mapped to primary seqs.''')
    parser.add_argument("-m", "--minscore", type=float, default=0, help='''PAF alignment file of assoc seqs mapped to primary seqs.''')
    parser.add_argument("-k", "--scorecolumn", type=int, default=12, help='''Integer. 1-based index of score column. Default = 12. Typical PAFs have MAPQ in column 12.''')

    args = parser.parse_args()
    return args




##############################################################################
''' FUNCTIONS '''
##############################################################################
    
def sortbed(bed, assume_bed_sorted=False):
    if not assume_bed_sorted:
        bed.sort(key=lambda x: (x[0], int(x[1])) )
    return bed

def convertPAFtoBEDlike(paf, fh=None, minscore=0, scorecolumn=11, mode='write', sortstoredbedlist=True):
    '''
    paf         = list object from PAF file made like this: [e.strip().split() for e in fh.readlines()]
    fh          = output file handle; either open connection in w mode or sys.stdout.
    minscore    = minimum value in scorecolumn to consider returning
    scorecolumn = column to find score in. Default = 11. (Mapq for minimap2).
    mode        = [ write | store | both ] = whether to write out or store and return as list object.
    '''
    assert mode in ('write', 'store', 'both')
    bedlikeList = []
    for aln in paf:
        name = aln[0]+"%"+aln[5]
        seqlen = aln[1]
        start = aln[2]
        end = aln[3]
        length = int(end) - int(start)
        mapq = aln[scorecolumn]
        matches = aln[9]
        alnLen = aln[10]
        pident = float(int(matches))/float(int(alnLen))
        pdiver = 1 - pident
        pmatch = float(int(matches))/float(length)
        dv = '-1'
        for ftr in aln:
            if ftr.startswith('dv:f:'): ## Approximate per-base sequence divergence based on minimizers
                dv = ftr.lstrip('dv:f:')
            if ftr.startswith('de:f:'): ## Gap-compressed per-base sequence divergence ; takes precedence over dv
                dv = ftr.lstrip('de:f:')
                break ## only put break here, not above, on purpose (to ensure de used over dv if both present)
        ambig = '0'
        for ftr in aln:
            if ftr.startswith('nn:i:'): ## Approximate per-base sequence divergence based on minimizers
                ambig = ftr.lstrip('nn:i:')
        alnLenUpdated = str(int(alnLen) - int(ambig))
        lengthUpdated = str(int(length) - int(ambig))
        pident2 = float(int(matches))/float(int(alnLenUpdated))
        pdiver2 = 1 - pident2
        pmatch2 = float(int(matches))/float(lengthUpdated)
        
        cv = 1-float(dv) ## convergence = complement of divergence
        
        ## OUT
        if float(mapq) >= minscore:
            bedline = [  str(  name ),
                         int(  start),
                         int(  end  ),
                         float(mapq ),
                         float(matches),
                         float(length),
                         float(alnLen),
                         float(pident),
                         float(pdiver),
                         float(pmatch),
                         float(ambig),
                         float(lengthUpdated),
                         float(alnLenUpdated),
                         float(pident2),
                         float(pdiver2),
                         float(pmatch2),
                         float(dv),
                         float(cv),
                         float(seqlen) ]
            if mode == 'write' or mode == 'both':
                fh.write( '\t'.join( str(out_e) for out_e in bedline ) + '\n' )
            if mode == 'store' or mode == 'both':
                bedlikeList.append( bedline )
    if mode == 'write':
        return None
    elif mode == 'store' or mode == 'both':
        if sortstoredbedlist:
            return sortbed( bedlikeList )
        else:
            return bedlikeList
        

def main():
    ## Parse args
    args = parse_args()
    scorecolumn = args.scorecolumn-1
    
    ### Get PAF
    with open(args.paf) as fh:
        paf = [e.strip().split() for e in fh.readlines()]

    ### CONVERT PAF TO BED-LIKE FMT
    fh = sys.stdout
    convertPAFtoBEDlike(paf = paf,
                        fh = fh,
                        minscore = args.minscore,
                        scorecolumn = scorecolumn)

##############################################################################
''' EXECUTE '''
##############################################################################

if __name__ == "__main__":
    main()
