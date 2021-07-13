#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="""

Given a fastq file - generate fasta and/or qual file

    """, formatter_class= argparse.RawTextHelpFormatter)



parser.add_argument('--fastq', '-f',
                   type= str, default='',
                   help='''Path to fastq file. For stdin, use stdin or -, or leave blank. ''')

parser.add_argument('--fa', action='store_true',
                   help='''Generate fasta file.''')

parser.add_argument('--qual', action='store_true',
                   help='''Generate qual file.''')

parser.add_argument('--out', '-o', type=str, default=False,
                   help='''Output prefix. Default: same as fastq. Use "stdout" or "-" to print to screen.''')

args = parser.parse_args()

assert args.fa or args.qual

if not args.out:
    args.out = args.fastq.split("/")[-1].split(".")[0]

if args.fastq in ('', '-', 'stdin'):
    args.fastq = sys.stdin

if args.out and args.out in ("stdout", "-"):
    out = sys.stdout
else:
    out = args.out+".fasta" if args.fa else args.out+".qual"

if args.fa:
    SeqIO.convert(args.fastq, "fastq", out, "fasta")
if args.qual:
    SeqIO.convert(args.fastq, "fastq", out, "qual")
