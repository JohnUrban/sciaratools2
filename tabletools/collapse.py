#!/usr/bin/env python3

import sys
import argparse
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser(description="""

Take in a table/txt file and collapse unique elements of 1 column
while doing operation on elements of a 2nd column.

No need for it to be sorted.
It will treat any occurence of an element in column 1 as part of same set.

    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--input', '-i',
                   type=str, required=True,
                   help='''Path to input text file''')

parser.add_argument('--delimiter', '-d',
                   type=str, default="\t",
                   help='''Delimiter. Default is tab.''')

parser.add_argument('--column','-c', type=int, default=1,
                    help=''' Column to collapse into unique instances..''')

parser.add_argument('--column2', '-c2', type=int, default=2,
                    help='''Column to perform collapsin operation on...''')

parser.add_argument('--operation', '-o', type=str, default='sum',
                    help=''' Operation to perform on column2. Default: sum.
Options: max, min, mean, median, sum, list, distinct.''')

parser.add_argument('--addcol', '-c3', type=str, default=False,
                    help='''Legacy (see --othercols too): When using ONLY max or min, also report these columns -- provide comma-separated list.''')

parser.add_argument('--othercols', '-c4', type=str, default=False,
                    help='''Collapse these other columns (comma-sep list) into their own columns in output as well. To get column of line for min or max see --addcol.''')

parser.add_argument('--otherOperation', '-O', type=str, default='list',
                    help='''Perform these operations (Options: max, min, mean, median, sum, list, distinct) on other columns from --othercols. If one provided, applied to all. Else need to provide same number of functions in comma-separated list.''')

parser.add_argument('--skip', '-s', type=str, default='#',
                    help='''Skip lines that start with given string. ''')

parser.add_argument('--header','-H', action='store_true', default=False,
                    help='''Print header line''')

parser.add_argument('--strings', '-str', action='store_true', default=False,
                    help=''' Default treats elements in column2 as floats. Treat "column2" as strings instead. Note: --operation must be set to list.... Doesn't apply to ADDCOLS and OTHERCOLS. This over-rides --mainType which defaults to float. However, this flag is somewhat redundant now since you could specify "--m str".''')

parser.add_argument('--mainType', '-m', type=str, default=float,
                    help='''How to interpret main column to collapse (float, int, str, etc). Default: float.''')

parser.add_argument('--otherType', '-M', type=str, default='str',
                    help='''How to interpret other columns to collapse (float, int, str, etc). Default: str. Applies to OTHERCOLS. If one provided, applied to all. Else need to provide same number of functions in comma-separated list.''')

args = parser.parse_args()

########################################################################################################
## FXNS
listfxn = lambda x: (",").join([str(e) for e in x])
distinctfxn = lambda x: (",").join([str(e) for e in list(set(x))])

########################################################################################################
## EXECUTE
if args.strings:
    assert args.operation in ["list", "distinct"]


########################################################################################################
## PYTHON INDEX OF COLUMNS 1 AND 2 (vols to collapse on and operate on)
c1 = args.column-1
c2 = args.column2-1

########################################################################################################
## GRAB FUNCTION(s) LIST
fxns = args.operation.split(",")


########################################################################################################
## DEFAULT ADD AND OTHER TO FALSE
addcol=False
othercols=False

########################################################################################################
## ADD COL ONLY WORKS WHEN APPLYING SAME MAIN FXN TO ALL ADDCOLS (AND FXN IS JUST MIN/MAX/LIST); ((some of this is just fossilized into the code from earlier versions, and I do recognize how hacky all this is.))
if args.addcol and len(fxns) == 1 and fxns[0] in ('min','max','list'):
    addcol=True
    ## "cols" with a "c" for ADDCOLS
    cols = [int(e)-1 for e in args.addcol.strip().split(",")]
    ## "D" is a dict with ADDCOL keys AND list vals.
    D = {}
    for col in cols:
        D[col] = defaultdict(list)

########################################################################################################
## OTHER COL IS ACTUALLY FOR WHEN YOU WANT TO SPECIFY OTHER FXNS FOR THE "ADDCOLS" ((so technically the main operation never need exceed 1))
##  THIS GETS PYINDEXES FOR OTHERCOLS; # OF OTHERCOLS; FXNS FOR OTHERCOLS; # FXNS ; 
if args.othercols:
    othercols=True
    ## "kols" with a "k" for OTHERCOLS
    kols = [int(e)-1 for e in args.othercols.strip().split(",")]
    nOtherKols = len(kols)

    ## Fxns for each othercol
    otherKolFxnTerms = args.otherOperation.strip().split(",")
    nOtherKolFxns = len(otherKolFxnTerms)
    assert nOtherKolFxns == 1 or nOtherKolFxns == nOtherKols

    ## Types for each othercol
    otherKolTypeTerms = args.otherType.strip().split(",")
    nOtherKolTypes = len(otherKolTypeTerms)
    assert nOtherKolTypes == 1 or nOtherKolTypes == nOtherKols

    ## Update term lists
    if nOtherKolFxns == 1:
        otherKolFxnTerms  = otherKolFxnTerms  * nOtherKols
    if nOtherKolTypes == 1:
        otherKolTypeTerms = otherKolTypeTerms * nOtherKols
     
    ## update lengths and assertions
    nOtherKolFxns  = len(otherKolFxnTerms)
    nOtherKolTypes = len(otherKolTypeTerms)
    assert nOtherKolFxns == nOtherKols


    ## gather functions
    otherKolFxns = []
    for e in otherKolFxnTerms:
        if e == 'sum':
            fxn = sum
        elif e == 'max':
            fxn = max
        elif e == 'min':
            fxn = min
        elif e == 'mean':
            fxn = np.mean
        elif e == 'median':
            fxn = np.median
        elif e == 'list':
            fxn = listfxn
        elif e == 'distinct':
            fxn = distinctfxn
        otherKolFxns.append(fxn)

    ## gather types
    otherKolTypes = []
    for e in otherKolTypeTerms:
        if e == 'str':
            tp = str
        elif e == 'float':
            tp = float
        elif e == 'int':
            tp = int
        elif e == 'list':
            tp = list
        elif e == 'tuple':
            fxn = tuple
        else:
            #print('DEBUG',e)
            assert e in ('str','int','float','list','tuple')
        otherKolTypes.append(tp)
        
        
    ## DD is a dict with OTHERCOL keys AND list vals; Each set of lines collapsed on MAINCOL will add to this (reset to empty for each set as part of processes below).
    DD = {}
    for col in kols:
        DD[col] = defaultdict(list)
##    if fxns[0] == 'max':
##        addfxn = max
##    elif fxns[0] == 'min':
##        addfxn = min


########################################################################################################
## GET OPERATION(S) FOR MAINCOL (AND ADDCOLS)
########################################################################################################
ops = []
for e in fxns:
    if e == 'sum':
        fxn = sum
    elif e == 'max':
        fxn = max
    elif e == 'min':
        fxn = min
    elif e == 'mean':
        fxn = np.mean
    elif e == 'median':
        fxn = np.median
    elif e == 'list':
        fxn = listfxn
    elif e == 'distinct':
        fxn = distinctfxn
    ops.append(fxn)


########################################################################################################
## CREATE DICT FOR INPUT PROCESSING
########################################################################################################
## "d" is a dict with MAINCOL keys AND list vals; 
d = defaultdict(list)


########################################################################################################
## ALLOW FOR STDIN OR FILE-BASED INPUT
########################################################################################################

if args.input == "-" or args.input == "stdin":
    f = sys.stdin
else:
    f = open(args.input)

########################################################################################################
## PROCESS INPUT LINES
########################################################################################################
for line in f:
    ## OPT. SKIP HEADER LINES
    if args.skip:
        if line.startswith(args.skip):
            continue

    ## PARSE LINE
    line = line.strip().split(args.delimiter)

    ## INTERPRETING MAIN_COL2 TYPES (MAIN OPERATION COL)... yes I realize the wonkiness of this all... its a monster built over many many years, with ornaments added every once in a while, sometimes redundantly in an effort to add a bell or whistle w/o wrapping head around what's already here!
    if args.strings:
        d[line[c1]].append(str(line[c2]))
    else:
        #d[line[c1]].append(float(line[c2]))
        d[line[c1]].append(args.mainType(line[c2]))

    ## 
    if addcol:
        for col in cols:
            D[col][line[c1]].append(line[col])
    if othercols:
        #for col in kols:
        #    DD[col][line[c1]].append( line[col] )
        ## Loop over indexes of # of OTHERCOLS
        for oi in range(len(kols)):
            ## Get the pyindex of current othercol
            col = kols[oi]

            ## coltype
            coltype = otherKolTypes[oi]

            ## interpreted
            interpreted_val = coltype(line[col])
            
            ## Interpret value and add
            DD[col][line[c1]].append( interpreted_val )



########################################################################################################
## PROCESS OUTPUT LINES
########################################################################################################

## OPT HEADER OUTPUT
if args.header:
    print((args.delimiter).join(["#c1_element", "number_found"] + fxns))


## LOOP OVER INPUT LINE MAINCOL VALUES
for e in sorted(d.keys()):
    ## Reset "ans" to empty list
    ans = []

    ## Append collapsed column2 (can do more than one op on col2; e.g. mean, min, max, list).
    for op in ops:
        ans.append( str(op(d[e])) )

    ## Note the number of lines that corresponded to this input value.
    length = len(d[e])

    ## ADDCOL ops ------------- this code is essentially deprecated; don't use --addcols and don't try to interpret any longer.
    if addcol:
        ## Get the index of the first value in ans....
        idx = d[e].index(float(ans[0]))
        ## reset addcols to empty list.
        addcols = []
        ## iterate over ADDCOLS
        for col in cols:
            ## Append __ to "addcols" list.
            addcols.append( D[col][e][idx] )
            
        print((args.delimiter).join([ str(e), str(length) ] + ans + addcols ))                   


    ## OTHERCOL ops on CURRENT INPUT LINE SET -- mutually exclusive from ADDCOL
    elif othercols:
        ## reset "othercols" to empty list.
        othercols = []
        #for col in kols:
        #    othercols.append( listfxn(DD[col][e]) )

        ## Loop over indexes of # of OTHERCOLS
        for oi in range(len(kols)):
            ## Get the pyindex of current othercol
            col = kols[oi]

            ## Get the function to use on current othercol
            applyfxn = otherKolFxns[oi]

            ##
            #print("DEBUG",col, e, applyfxn)
            #print("DEBUG",DD[col][e])
            newres = applyfxn(DD[col][e])
            

            ## Append new result to "othercols" list.
            othercols.append( str(newres) )
        
        print((args.delimiter).join([ str(e), str(length) ] + ans + othercols ))                   

    else:    
        print((args.delimiter).join([ str(e), str(length) ] + ans))
        

















