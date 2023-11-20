#!/usr/bin/env python3

import sys, argparse
from collections import defaultdict



def parse_args():
    parser = argparse.ArgumentParser(description="""

    DESCRIPTION -

        Given table,
            keep only 1 line for each element in name column 
            that has the highest score in score column.
            
            

        Use case:
        - BLASTP output on gene set, keep entry with highest bit score.
        - BED interval with highest (or lowest) end coordinate per chromosome.
        

        """, formatter_class= argparse.RawTextHelpFormatter)

    parser.add_argument('table', metavar='table', nargs='*',
                       type= str, 
                       help='''Paths to as many tables as you want. They will be treated as one giant table.
                            This can also be left empty for standard in or specified as - or stdin.''')

    parser.add_argument('-n', '--name_column', type=int, default=1,
                        help='''Column where elements/names are found.''')

    parser.add_argument('-s', '--score_column', type=int, default=2,
                        help='''Column where scores are found.''')


    parser.add_argument('-d', '--delim', type=str, default='\t',
                        help='''Delimiter. Default = tab.''')

    parser.add_argument('-L', '--lowest', action='store_true', default=False,
                        help='''Filter for lowest value instead of highest.''')

    args = parser.parse_args()
    return args


def tableFilter_step01(tables):
    ''' tables is a list of filenames for tabular-like files'''
    if len(tables) == 0 or tables[0] in ('stdin', '-') or tables[0].startswith('<('):
        tables = [sys.stdin]
        def opentable(x):
            return x
        def closetable(x):
            return x
    else:
        def opentable(x):
            return open(x)
        def closetable(x):
            return x.close()
    return tables, opentable, closetable


def tableFilter_step02(tables, delim, score_col, name_col, opentable, closetable, lowest=False):
    ''' tables is a list of filenames for tabular-like files'''
    lines = {}
    scores = {}
    order = []
    for table in tables:
        tablefile = opentable(table)
        for line in tablefile:
            line = line.strip().split(delim)
            score = float(line[score_col])
            order.append(line[name_col])
            try:
                if score > scores[line[name_col]] and not lowest:
                    scores[line[name_col]] = score
                    lines[line[name_col]] = line
                elif score < scores[line[name_col]] and lowest:
                    scores[line[name_col]] = score
                    lines[line[name_col]] = line
            except:
                scores[line[name_col]] = score
                lines[line[name_col]] = line
        closetable(tablefile)
    return lines, scores, order
    

def tableFilter_step03(lines, order, delim):
    ''' tables is a list of filenames for tabular-like files'''
    ## Output in order of appearance in input
    for name in order:
        try:
            print( (delim).join(lines[name]) )
            lines.pop(name)
        except:
            ## The intention is to cause it to fail silently after seeing element once and printing
            ## However, this can cause it to fail silently in unanticipated scenarios as well.
            pass

def tableFilterPipeline(tables, delim, score_col, name_col, lowest):
    ''' tables is a list of filenames for tabular-like files'''
    ## Step 1
    tables, opentable, closetable = tableFilter_step01( tables = tables )

    ## Step 2
    lines, scores, order = tableFilter_step02(  tables      = tables,
                                                delim       = delim,
                                                score_col   = score_col,
                                                name_col    = name_col,
                                                opentable   = opentable,
                                                closetable  = closetable,
                                                lowest      = lowest)
    ## Step 3
    tableFilter_step03(lines, order, delim)

    ## Return None
    return None
    

def main():
    ## Get ARGs
    args = parse_args()

    ## Process ARGs
    score_col = args.score_column - 1
    name_col = args.name_column - 1

    ## Execute
    tableFilterPipeline(tables      = args.table,
                        delim       = args.delim,
                        score_col   = score_col,
                        name_col    = name_col,
                        lowest      = args.lowest)


    


##############################################################################
''' EXECUTE '''
##############################################################################

if __name__ == "__main__":
    main()







##############################################################################
''' DEPRECATED CODE '''
##############################################################################


## BELOW DOES NOT WORK WITH STDIN
##for table in args.table:
##    tablefile = opentable(table)
##    for line in tablefile:
##        line = line.strip().split(args.delim)
##        try:
##            print (args.delim).join(lines[line[name_col]])
##            lines.pop(line[name_col])
##        except:
##            ## The intention is to cause it to fail silently after seeing element once and printing
##            ## However, this can cause it to fail silently in unanticipated scenarios as well.
##            pass
##    closetable(tablefile)


