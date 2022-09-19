#!/usr/bin/env python3

import sys, argparse
from Bio import SeqIO



def parse_args():

    parser = argparse.ArgumentParser(description="""

    Given a fasta file, return sequences with ends stripped according to instructions.

    Main use case:
    - Remove flanking NNNNNNN sequences from contigs and scaffolds.

    Other use cases:
    - very basic trimming of reads
        - recommended to use tools designed for this application instead, e.g.:
            - fastp
            - trimmomatic

        """, formatter_class= argparse.RawTextHelpFormatter)

    parser.add_argument('fastx', metavar='fastx', nargs='?',
                        type= str, default="-",
                        help='''Path to fasta/fastq. For stdin, leave blank or use 'stdin' or '-', or the '<(cmds) approach. --fa/--fq must be specified with stdin.''')
    ## ? = 0 or 1 pos args ; if 0, default
    ## * = 0 or more in list
    ## + = 1 or more in list ; 0 throws error

    filetype = parser.add_mutually_exclusive_group()
    filetype.add_argument('--fa', action='store_true',
                          help='''Explicitly define as fasta input. --fa/--fq must be specified with stdin.''')
    filetype.add_argument('--fq', action='store_true',
                          help='''Explicitly define as fastq input. --fa/--fq must be specified with stdin.''')

    parser.add_argument("-c", "--flankchar", type=str, default='N', help='''Character or pattern to strip from flanks. Default = N. See --ignorecase (e.g. for N or n).''')
    parser.add_argument("-i", "--ignorecase", action='store_true', default=False, help='''Try harder to strip flankchar or pattern by also stripping char.upper() and char.lower()''')
    
    args = parser.parse_args()

    ## PROCESS ARGS
    is_stdin    = args.fastx in ('stdin','-') or args.fastx.startswith('<(')
    ftype_known = args.fa or args.fq
    if is_stdin and not ftype_known:
        print("When using stdin, define filetype with --fa or --fq.")
        quit()


    ## if fa/fq not specified, determine:
    if not ftype_known:
        is_fq = args.fastx.endswith('fq') or args.fastx.endswith('fastq')
        is_fa = args.fastx.endswith('fa') or args.fastx.endswith('fasta') or args.fastx.endswith('fna') or args.fastx.endswith('fpa')
        if is_fq:
            args.fq = True
        elif is_fa:
            args.fa = True
        else:
            with open(args.fastx) as fh:
                line = next(fh)
                if line.startswith('>'):
                    args.fa = True
                elif line.startswith('@'):
                    args.fq = True
                else:
                    print("Unable to deterimine filetype. Specify as --fa/--fq.")
                    quit()
    if args.fa:
        args.x = "fasta"
    elif args.fq:
        args.x = "fastq"

    if is_stdin:
        args.fastx = sys.stdin

    
    return args

def strip(record, flank='N', ignorecase=False):
    if ignorecase:
        return ">"+ record.description + "\n" + str(record.seq).strip().strip(flank).strip(flank.upper()).strip(flank.lower()) + "\n"
    else:
        return ">"+ record.description + "\n" + str(record.seq).strip().strip(flank) + "\n"



def main():
    args = parse_args()
    for record in SeqIO.parse(args.fastx, args.x):
        sys.stdout.write( strip(record      = record,
                                flank       = args.flankchar,
                                ignorecase  = args.ignorecase) )










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

