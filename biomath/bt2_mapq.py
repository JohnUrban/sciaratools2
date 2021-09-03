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
    This has mainly been tested for single-end reads in end-to-end mode at the moment...
    When written at least, bowtie2 assigned scores between 0-42.
    I now also see MAPQ=44 in paired-end datasets....
    ....possibly when both mates would have received 42.

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
parser.add_argument("-XS", "--XS", type=str, default=None, help='''Alignment score of next best alignment -- needs to be the same or lower than AS. Do not use if no XS given. Or use "-1" or "None".''')
parser.add_argument("-f", "--table", type=str, default=None, help='''Path to table file with AS and XS columns''')
parser.add_argument("-ASC", "--AScolumn", type=int, default=None, help='''1-based Column to find AS in if table file given.''')
parser.add_argument("-XSC", "--XScolumn", type=int, default=None, help='''1-based Column to find XS in if table file given.''')
parser.add_argument("-ASC2", "--AScolumn2", type=int, default=None, help='''Tables that describe a pair per line. 1-based Column to find AS of mate in if table file given.''')
parser.add_argument("-XSC2", "--XScolumn2", type=int, default=None, help='''Tables that describe a pair per line. 1-based Column to find XS of mate in if table file given.''')
parser.add_argument('-r', '--readlength', type=int, default=None, help=''' Scores are often dependent on read length. Thus, it must be specified here.''')
parser.add_argument("-s", "--scoremin", type=str, default=None, help=''' Argument to bowtie2 --score-min. Default: None. Use --endtoend or --local to specify defaults of "L,-0.6,-0.6" or "G,20,8" respectively. See bowtie2 documentation for more info: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml''')
parser.add_argument("-ma", "--matchscore", type=float, default=None, help='''"--ma" tells bowtie2 how to score a match. The defaults in endtoend and local modes are 0 and 2, respectively. If defaults not used in bowtie2, this needs to be changed to what --ma was set to! This is used to calculate the perfect score given readlength. ''')
default = parser.add_mutually_exclusive_group(required=True)
default.add_argument("-e", "--endtoend", action='store_true', default=None, help='''EndtoEnd and Local scores are calculated differently, so this needs to be specified. Also helps specify default --score-min if bowtie2 was run in --end-to-end mode w/o changing the --score-min default therein: L,-0.6,-0.6. The --scoremin option overrides this.''')
default.add_argument("-l", "--local", action='store_true', default=None, help='''EndtoEnd and Local scores are calculated differently, so this needs to be specified. Also helps specify default --score-min if bowtie2 was run in --local mode w/o changing the --score-min default therein: G,20,8. The --scoremin option overrides this.''')
args = parser.parse_args()



################################################################
''' FUNCTIONS '''
################################################################

def endtoend_MAPQ(AS, XS, scMin, diff, bestOver, bestdiff):
    ## If does not have second best
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
    else:   ## it does have second best
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

            
def local_MAPQ(AS, XS, scMin, diff, bestOver, bestdiff):
    ## BELOW - I have changed ">=" in original bowtie2 code to ">" where it gave more results consistent w/ bt2 output
    ##      Future testing may lead me to change it back, but it fixed hundreds of discrepancies in original testing of paired MAPQ calculations
    ## If does not have second best
    if XS < scMin: ## if only single alignment (AS and no XS)
        if bestOver > diff*0.8:
            return 44 # 44 is best in local mode, but only 42 in endtoend mode...
        elif bestOver >= diff*0.7:
            return 42
        elif bestOver > diff*0.6: 
            return 41
        elif bestOver >= diff*0.5:
            return 36
        elif bestOver >= diff*0.4: 
            return 28
        elif bestOver >= diff*0.3:
            return 24
        else:
            return 22
    else:   ## it does have second best
        if bestdiff >= diff*0.9:
            return 40
        elif bestdiff > diff*0.8:
            return 39
        elif bestdiff >= diff*0.7:
            return 38
        elif bestdiff > diff*0.6: 
            return 37
        elif bestdiff >= diff*0.5:
            if bestOver == diff:
                return 35
            elif bestOver >= diff*0.5:
                return 25
            else:
                return 20
        elif bestdiff > diff*0.4: 
            if bestOver == diff:
                return 34
            elif bestOver >= diff*0.5:
                return 21
            else:
                return 19
        elif bestdiff > diff*0.3:
            if bestOver == diff:
                return 33
            elif bestOver >= diff*0.5:
                return 18
            else:
                return 16
        elif bestdiff > diff*0.2:
            if bestOver == diff:
                return 32
            elif bestOver >= diff*0.5:
                return 17
            else:
                return 12
        elif bestdiff > diff*0.1:
            if bestOver == diff:
                return 31
            elif bestOver >= diff*0.5:
                return 14
            else:
                return 9
        elif bestdiff > 0:
            if bestOver >= diff*0.5:
                return 11
            else:
                return 2
        else: ## best and 2nd best AS are the same (multireads where AS=XS)
            if bestOver >= diff*0.5: 
                return 1
            else:
                return 0
    
def bt2_mapq(AS, XS=None, alnmode=None, scMin=None, scPer=None):
    '''
    scMin = minScore
    scPer = perfect score
    '''
    
    ## Did the read have a second-best alignment?
    if XS == None: ## make it less than scMin
        ## scMin = score of a just barely valid match
        XS = scMin-1
    #if XS > AS:
    #    return None

    ## Difference between the perfect score and the minimum score
    ## range of aln scores and of largest possible diff between XS-AS
    diff = max(1, abs(scPer-scMin))

    ## Best alignment score found
    best = max(AS, XS)
    
    ## Difference between best alignment score seen and score minumum
    ## bestOver = AS-scMin
    bestOver = best-scMin ## Gives diff between AS and minScore as a positive
                        ## this is biggest that bestdiff (AS-XS) can be -- after AS is seen


    ## Absolute difference between the best and second best alignment scores (usually AS>XS)
    if AS > XS:
        bestdiff = abs(abs(AS)-abs(XS)) ## gives diff in AS between best and 2nd best -- biggest diff after both AS and XS seen
    else:
        bestdiff = 1 ## seems like 0 or negative would be better, but so far this has minimized discrepancies

    
    ## Was alignment mode in end-to-end or local?
    if alnmode == 'endtoend':
        return endtoend_MAPQ(AS, XS, scMin, diff, bestOver, bestdiff)
    elif alnmode == 'local':
        return local_MAPQ(AS, XS, scMin, diff, bestOver, bestdiff)
    else:
        ## An error catcher here, but the error should be impossible
        sys.stderr.write('MAPQ calculation error: alnmode must be endtoend or local.\n')
        quit()
    
    



def getscoremin(fxn, readlength):
    l = fxn.strip().split(',')
    fd = {'C':lambda b,m,x: b + m,
          'G':lambda b,m,x: b + m*np.log(x),
          'L':lambda b,m,x: b + m*x,
          'S':lambda b,m,x: b + m*np.sqrt(x)}
    return fd[l[0]](float(l[1]), float(l[2]), int(readlength))

def run_simple(args, alnmode, scoremin, scoreperfect):
    args.XS = None if args.XS in (None, "None", -1) else int(args.XS)
    mapq = bt2_mapq(
                    AS=args.AS,
                    XS=args.XS,
                    alnmode=alnmode,
                    scMin=scoremin,
                    scPer=scoreperfect)
    out = '\t'.join(['AS', 'XS', 'ScoreMin', 'MAPQ']) + '\n'
    out += '\t'.join(str(e) for e in [args.AS,
                                      args.XS,
                                      scoremin,
                                      mapq])  + '\n'
    sys.stdout.write( out )


def run_table(args, alnmode, scoremin, scoreperfect):
    asc = args.AScolumn - 1
    xsc = args.XScolumn - 1
    paired=False
    if args.AScolumn2 is not None and args.XScolumn2 is not None:
        paired = True
        asc2 = args.AScolumn2 - 1
        xsc2 = args.XScolumn2 - 1
    #Get ncols
    with open(args.table) as tmp:
        line = next(tmp).strip().split()
        ncol = len(line)
        
    if not paired:
        out = '\t'.join(['ucol_'+str(i) for i in range(ncol)] + ['MAPQ']) + '\n'
    else:
        out = '\t'.join(['ucol_'+str(i) for i in range(ncol)] + ['MAPQ1', 'MAPQ2', 'J1', 'J2', 'J3', 'J4', 'J5', 'J6']) + '\n'

    with open(args.table) as table:
        for row in table:
            row = row.strip().split()
            AS = float(row[asc])
            XS = None if row[xsc] in (None, "None", -1, "-1") else float(row[xsc])
            row += ["MAPQ"]
            mapq1 = bt2_mapq(
                    AS=AS,
                    XS=XS,
                    alnmode=alnmode,
                    scMin=scoremin,
                    scPer=scoreperfect)
            row += [mapq1]
            if paired:
                AS2 = float(row[asc2])
                XS2 = None if row[xsc2] in (None, "None", -1, "-1") else float(row[xsc2])
                mapq2 = bt2_mapq(
                    AS=AS2,
                    XS=XS2,
                    alnmode=alnmode,
                    scMin=scoremin,
                    scPer=scoreperfect)
                row += [mapq2]
                
                row += ["JOINT"]

                ## JOINT1
                ASJ = AS + AS2

                ## all below give same result; combining 1 with if statements does a little better...
                #1 XSJ = None if (XS is None and XS2 is None) else ((XS if XS is not None else 0) + (XS2 if XS2 is not None else 0))
                #2 XSJ = None if (XS is None and XS2 is None) else ((XS if XS is not None else scoremin-1) + (XS2 if XS2 is not None else scoremin-1))
                #3 XSJ = (scoremin-1)*2 if (XS is None and XS2 is None) else ((XS if XS is not None else scoremin-1) + (XS2 if XS2 is not None else scoremin-1))
                #4 XSJ = (scoremin-1)*2 if (XS is None and XS2 is None) else ((XS if XS is not None else scoremin) + (XS2 if XS2 is not None else scoremin))


                ## combining 1 with if statements does a little better
                #XSJ = None if (XS is None and XS2 is None) else ((XS if XS is not None else 0) + (XS2 if XS2 is not None else 0))
                #if XSJ is not None:
                #    if XSJ == XS:
                #        XSJ = XS + max(AS, AS2) 
                #    elif XSJ == XS2:
                #        XSJ = XS2 + max(AS, AS2) 
                #if XSJ is None:
                #    XSJ = (scoremin-1) * 2

                ##Below is similar to the above code, but the more-sane human readable syntax allowed me to optimize a bit further
                if XS is not None and XS2 is None:
                    ##XSJ = XS + max(AS, AS2)
                    XSJ = XS + AS2
                elif XS is None and XS2 is not None:
                    ##XSJ = XS2 + max(AS, AS2)
                    XSJ = AS + XS2
                elif XS is None and XS2 is None:
                    XSJ = (scoremin-1) * 2 ## This can simply be scoremin and still work....at least in localmode C,40,0
                    # XSJ = min(AS, AS2) * 2            ....no
                    # XSJ = min(AS, AS2) + scoremin-1    ....no
                else: ## neither none
                    XSJ = XS + XS2
                    ## XSJ = max(XS, XS2)*2 ; when the above fails to give correct pair mapq, this at least sometimes solves it

                # The pair XS sometimes exceeds the pair of AS, but very rarely, and correcting for it did not further optimize
                #if XSJ > ASJ:
                #    XSJ = ASJ


                mapqJ1 = bt2_mapq(
                    AS=ASJ,
                    XS=XSJ,
                    alnmode=alnmode,
                    scMin = scoremin * 2,
                    scPer = scoreperfect * 2
                    )
                row += [mapqJ1]

##                ## JOINT2
##                ASJ = 0.5*sum([AS,AS2])
##                XSJ = None if (XS is None and XS2 is None) else 0.5*((XS if XS is not None else 0) + (XS2 if XS2 is not None else 0))
##                mapqJ2 = bt2_mapq(
##                    AS=ASJ,
##                    XS=XSJ,
##                    alnmode=alnmode,
##                    scMin = scoremin,
##                    scPer = scoreperfect)
##                row += [mapqJ2]
##
##                ## JOINT3
##                ASJ = np.sqrt(AS * AS)
##                XSJ = None if (XS is None and XS2 is None) else np.sqrt((XS if XS is not None else 1) * (XS2 if XS2 is not None else 1))
##                mapqJ3 = bt2_mapq(
##                    AS=ASJ,
##                    XS=XSJ,
##                    alnmode=alnmode,
##                    scMin = scoremin,
##                    scPer = scoreperfect)
##                row += [mapqJ3]
##
##                ## JOINT4
##                mapqJ4 = round(0.5*((0 if mapq1 is None else mapq1) + (0 if mapq2 is None else mapq2)))
##                row += [mapqJ4]
##
##                ## JOINT5
##                mapqJ5 = round(np.sqrt((1 if mapq1 is None else mapq1)*(1 if mapq2 is None else mapq2)))
##                row += [mapqJ5]
##
##                ## JOINT6
##                ASJ = AS + AS2
##                XSJ = None if (XS is None and XS2 is None) else ((XS if XS is not None else 0) + (XS2 if XS2 is not None else 0))
##                mapqJ1 = bt2_mapq(
##                    AS=ASJ,
##                    XS=XSJ,
##                    alnmode=alnmode,
##                    scMin = scoremin,
##                    scPer = scoreperfect)
##                row += [mapqJ1]
##
##                ## JOINT5
##                mapqJ5 = round(np.log(
##                    np.exp((0 if mapq1 is None else mapq1) + (0 if mapq2 is None else mapq2))))
##                row += [mapqJ5]

            out += '\t'.join(str(e) for e in row) + '\n'
    sys.stdout.write( out )

def run(args):
    simple_exe = (args.AS is not None)
    table_exe = (args.table is not None and args.AScolumn is not None and args.XScolumn is not None)
    scoreminfxnprovided = (args.scoremin is not None) or (args.local is not None) or (args.endtoend is not None)
    readlengthgiven = args.readlength is not None
    assert not (simple_exe and table_exe)
    assert (simple_exe or table_exe)
    assert scoreminfxnprovided
    assert readlengthgiven


    ## Get aln mode
    if args.endtoend:
        alnmode = "endtoend"
    elif args.local:
        alnmode = "local"

    ## Getscoremin default function if None provided.
    if args.scoremin is None:
        if args.local:
            args.scoremin = 'G,20,8'
        elif args.endtoend:
            args.scoremin = 'L,-0.6,-0.6'

    ## Get minimum valid aln score (--score-min) given scoremin function/parameters, and read length
    scoremin = getscoremin(fxn = args.scoremin,
                           readlength = args.readlength)

    ## Get maximum valid aln score -- i.e. the perfect score.
    if args.matchscore is None:
        if args.endtoend:
            args.matchscore = 0.0
        elif args.local:
            args.matchscore = 2.0
    scoreperfect = args.matchscore * args.readlength


    
    if simple_exe:
        run_simple(args, alnmode, scoremin, scoreperfect)
    elif table_exe:
        run_table(args, alnmode, scoremin, scoreperfect)
    else:
        pass ## this should never happen

    ## Return ans
    

################################################################
''' EXECUTE '''
################################################################
run(args)

