#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, sys
import numpy as np

## transcribed from C code at:
## http://seqanswers.com/forums/showpost.php?p=119928&postcount=510
## described in detail at:
##  http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html


# Parse command-line arguments
parser = argparse.ArgumentParser(description='''
Description:
    This really only works / been tested for end-to-end mode at the moment...

    Computes estimated MAPQ from bowtie2 given an alignment score (AS) and an optional second alignment score (XS) that is the same or lower than AS.

    Can also provide a table file, specifying the columns for each.

    For paired-end reads,
        - alignment scores are calculated as the sum of the AS for each mate

    
    Bowtie2 notes:

    
    - --score-min F,B,M :: F=transform function (C=constant, G=log, L=linear, S=sqrt); B=b in mx+b; M=m in mx+b
        - computed min_score(read_length) = B + M*F(read_length)
        - Minimum score = absolute minimum + length-dependent addition
        - If the "absolute minimum" is set to 0, then it is purely dependent on read length.
        - Otherwise, it requires first meeting a constant minimum like BWA, followed by becoming more stringent with read length.

    - Default for end-to-end mode = L,-0.6,-0.6
        - L,-0.6,-0.6(read_length) = -0.6 + -0.6*read_length
        - Read_length -> MinScore
            - 20 -> -12.6
            - 30 -> -18.6
            - 50 -> -30.6
            - 75 -> -45.6
            - 100 -> -60.6
            - 150 -> -90.6
            
    - Default for local mode = G,20,8
        - G,20,8(read_length) = 20 + 8*log(read_length)
        - Read_length -> MinScore
            - 20 -> 43.9659
            - 30 -> 47.2096
            - 50 -> 51.2962
            - 75 -> 54.5399
            - 100 -> 56.8414
            - 150 -> 60.0851
            
    - Purely dependent on read length:
        - End-to-end mode - similar to default without absolute minimum:
            - L,0,-0.6
        - Local mode - similar to default without absolute minimum:
            - G,0,8
            
    - Only requires passing absolute minimum / independent of read length
        - In BWA, the alignment score is calculated independent of read length, and the default minscore is 30 (-T 30).
        - If you wanted Bowtie2 to have a similar approach, the following codes work by analogy (Bowtie2 uses +2 for matches by default, BWA +1):
            - -T 20 ==> C,40,0
                - C,40,0(read_length) = 40 + 0*read_length = 40
            - -T 30 ==> C,60,0
                - C,60,0(read_length) = 60 + 0*read_length = 60
            - ETC
        - If you wanted it even more similar, the following should do the trick in theory, but differences are found anyway:
            - -T 20 ==> --ma 1 --mp 4,4 --rdg 6,1 --rfg 6,1 --np 4 --score-min C,20,0
            - -T 30 ==> "--ma 1 --mp 4,4 --rdg 6,1 --rfg 6,1 --np 4 --score-min C,30,0
            - ETC


    Code that transforms alignment scores to mapq:
    ## transcribed from C code at:
    ##      http://seqanswers.com/forums/showpost.php?p=119928&postcount=510
    ## described in detail at:
    ##      http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html



''', formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-AS", "--AS", type=int, default=None, help='''Alignment score.''')
parser.add_argument("-XS", "--XS", type=int, default=None, help='''Alignment score of next best alignment -- needs to be the same or lower than AS.''')
parser.add_argument("-f", "--table", type=str, default=None, help='''Path to table file with AS and XS columns''')
parser.add_argument("-ASC", "--AScolumn", type=int, default=None, help='''1-based Column to find AS in if table file given.''')
parser.add_argument("-XSC", "--XScolumn", type=int, default=None, help='''1-based Column to find XS in if table file given.''')
parser.add_argument('-r', '--readlength', type=int, default=None, help=''' Scores are often dependent on read length. Thus, it must be specified here.''')
parser.add_argument("-s", "--scoremin", type=float, default=None, help=''' Argument to bowtie2 --score-min. Default: None. Use --endtoend or --local to specify defaults of "L,-0.6,-0.6" or "G,20,8" respectively. See bowtie2 documentation for more info: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml''')
default = parser.add_mutually_exclusive_group()
default.add_argument("-e", "--endtoend", action='store_true', default=None, help='''Helps specify default --score-min if bowtie2 was run in --end-to-end mode w/o changing the --score-min default therein: L,-0.6,-0.6''')
default.add_argument("-l", "--local", action='store_true', default=None, help='''Helps specify default --score-min if bowtie2 was run in --local mode w/o changing the --score-min default therein: G,20,8''')
args = parser.parse_args()



################################################################
''' FUNCTIONS '''
################################################################
def bt2_mapq_end2end(AS, XS=None, scMin=-30.6):
    '''scMin = minScore'''
    if XS == None: ## make it less than scMin
        XS = scMin-1
    if XS > AS:
        return None
    diff = abs(scMin) ## range of aln scores and of largest possible diff between XS-AS before either is seen
    bestOver = AS-scMin ## Gives diff between AS and minScore as a positive
                        ## this is biggest that bestdiff (AS-XS) can be -- after AS is seen
    bestdiff = abs(abs(AS)-abs(XS)) ## gives diff in AS between best and 2nd best -- biggest diff after both AS and XS seen
    if XS < scMin: ## if only single alignment (no XS)
        if bestOver >= diff*0.8:
            return 42 #42 when 
        elif bestOver >= diff*0.7:
            return 40
        elif bestOver >= diff*0.61: ## originally >= 0.60 but this made it agree with my bt2 tests better
            return 24
        elif bestOver >= diff*0.5:
            return 23
        elif bestOver >= diff*0.42: ## originally >= 0.40 but this made it agree with my bt2 tests better
            return 8
        elif bestOver >= diff*0.3:
            return 3
        else:
            return 0
    else:
        if bestdiff >= diff*0.9:
            if bestOver == diff:
                return 39
            else:
                return 33
        elif bestdiff >= diff*0.8:
            if bestOver == diff:
                return 38
            else:
                return 27
        elif bestdiff >= diff*0.97:
            if bestOver == diff:
                return 37
            else:
                return 26
        elif bestdiff >= diff*0.6:
            if bestOver == diff:
                return 36
            else:
                return 22
        elif bestdiff >= diff*0.5:
            if bestOver == diff:
                return 35
            elif bestOver >= diff*0.84:
                return 25
            elif bestOver >= diff*0.68:
                return 16
            elif bestOver >= diff*0.68:
                return 5
        elif bestdiff >= diff*0.4:
            if bestOver == diff:
                return 34
            elif bestOver >= diff*0.84:
                return 21
            elif bestOver >= diff*0.68:
                return 14
            else:
                return 4
        elif bestdiff >= diff*0.3:
            if bestOver == diff:
                return 32
            elif bestOver >= diff*0.88:
                return 18
            elif bestOver >= diff*0.67:
                return 15
            else:
                return 3
        elif bestdiff >= diff*0.2:
            if bestOver == diff:
                return 31
            elif bestOver >= diff*0.88:
                return 17
            elif bestOver >= diff*0.67:
                return 11
            else:
                return 0
        elif bestdiff >= diff*0.1:
            if bestOver == diff:
                return 30
            elif bestOver >= diff*0.88:
                return 12
            elif bestOver >= diff*0.67:
                return 7
            else:
                return 0
        elif bestdiff > 0:
            if bestOver >= diff*0.67:
                return 6
            else:
                return 2
        else: ## best and 2nd best AS are the same (multireads where AS=XS)
            if bestOver >= diff*0.68: ## originally >= 0.67 but this made it agree with my bt2 tests better
                return 1
            else:
                return 0



def getscoremin(fxn, readlength):
    l = fxn.strip().split(',')
    fd = {'C':lambda b,m,x: b + m,
          'G':lambda b,m,x: b + m*np.log(x),
          'L':lambda b,m,x: b + m*x,
          'S':lambda b,m,x: b + m*np.sqrt(x)}
    return fd[l[0]](float(l[1]), float(l[2]), int(readlength))

def run_simple(args, scoremin):
    mapq = bt2_mapq_end2end(
                    AS=args.AS,
                    XS=args.XS,
                    scMin=scoremin)
    out = '\t'.join(['AS', 'XS', 'ScoreMin', 'MAPQ']) + '\n'
    out += '\t'.join(str(e) for e in [args.AS,
                                      args.XS,
                                      scoremin,
                                      mapq])  + '\n'
    sys.stdout.write( out )


def run_table(args, scoremin):
    asc = args.AScolumn - 1
    xsc = args.XScolumn - 1
    out = '\t'.join(['AS', 'XS', 'MAPQ']) + '\n'
    with open(args.table) as table:
        for row in table:
            row = row.strip().split()
            AS = float(row[asc])
            XS = float(row[xsc])
            mapq = bt2_mapq_end2end(
                    AS=AS,
                    XS=XS,
                    scMin=scoremin)
            row += [mapq]
            out += '\t'.join(str(e) for e in row) + '\n'
    sys.stdout.write( out )

def run(args):
    simple_exe = (args.AS is not None and args.XS is not None)
    table_exe = (args.table is not None and args.AScolumn is not None and args.XScolumn is not None)
    scoreminfxnprovided = (args.scoremin is not None) or (args.local is not None) or (args.endtoend is not None)
    readlengthgiven = args.readlength is not None
    assert not (simple_exe and table_exe)
    assert (simple_exe or table_exe)
    assert scoreminfxnprovided
    assert readlengthgiven


    ## Getscoremin function
    if args.scoremin is None:
        if args.local:
            args.scoremin = 'G,20,8'
        elif args.endtoend:
            args.scoremin = 'L,-0.6,-0.6'

    ## Get score min
    scoremin = getscoremin(fxn = args.scoremin,
                           readlength = args.readlength)
    
    if simple_exe:
        run_simple(args, scoremin)
    elif table_exe:
        run_table(args, scoremin)
    else:
        pass ## this should never happen

    ## Return ans
    

################################################################
''' EXECUTE '''
################################################################
run(args)

