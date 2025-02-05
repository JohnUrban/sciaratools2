#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from collections import Counter
#import pandas as pd

parser = argparse.ArgumentParser(description="""

DESCRIPTION

    THIS TOOL HAS NOT BEEN DEVELOPED AND SHOULD NOT BE USED.
    
    =======================
    - Take a multi-fasta file (.fa file with >=1 entry) where sequences may contain non-canonical or ambiguous amino acid letters.
    - Output report.

    
    Anticipated canonical letters are in "A,a,C,c,D,d,E,e,F,f,G,g,H,h,I,i,K,k,L,l,M,m,N,n,P,p,Q,q,R,r,S,s,T,t,V,v,W,w,Y,y"
    - AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYy

    Anticipated non-canonical letters are in: OoUu
        - these are converted to a random canonical AA 

    Anticipated ambiguous letters are in BbJjXxZz
        - B converted to D or N
        - J converted to I or L
        - X converted to any canonical AA
        - Z converted to E or Q
    

    Unanticipated letters are any other character encountered in sequence. These are treated like "X" by default.






    Canonical AA Table
    ===================
    3-ltr 1-ltr Name
    ------------------
    Ala	A Alanine
    Arg	R Arginine
    Asn	N Asparagine
    Asp	D Aspartic acid
    Cys	C Cysteine
    Gln	Q Glutamine
    Glu	E Glutamic acid
    Gly	G Glycine
    His	H Histidine
    Ile	I Isoleucine
    Leu	L Leucine
    Lys	K Lysine
    Met	M Methionine
    Phe	F Phenylalanine
    Pro	P Proline
    Ser	S Serine
    Thr	T Threonine
    Trp	W Tryptophan
    Tyr	Y Tyrosine
    Val	V Valine


    Non-canonical AA Table
    ===================
    3-ltr 1-ltr Name
    ------------------
    Pyl	O Pyrrolysine
    Sec	U Selenocysteine


    Ambiguous AA Table
    ===================
    3-ltr 1-ltr Name
    ------------------
    Asx	B Aspartic acid or Asparagine
    Glx	Z Glutamic acid or Glutamine
    Xaa	X Any amino acid
    Xle	J Leucine or Isoleucine
    
    """, formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument('--fasta', '-f', required=True,
                   type=str, default=False,
                   help='''Provide path to fasta file''')

parser.add_argument('--outputType', '-O', 
                   type=str, default='report',
                   help='''Options are 'report' (default) or 'table'.''')

args = parser.parse_args()



safe_ltrs = "AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYy"
other_ltrs = "BbJjXxZzOoUu"

## functions
def probs(ltrs, freq):
    s = 0
    for b in ltrs:
        s += freq[b]
    p = []
    for b in ltrs:
        p.append( freq[b]/s )
    return p

def report(fa, freq):
    print(fa.name)
    length = sum(freq.values())

    
    #for key in freq.keys():
    for key, v in sorted(freq.items(), key=lambda item: item[1], reverse=True):
        print(key + " " + str(round(100*v/length, 2))+"%")
    print("Seq Length: " + str(len(fa)))
    print("Amino Acids profiled: " + " " + str(length))
    print('========================================================')



for fa in SeqIO.parse(args.fasta, 'fasta'):
    freq = Counter()
    for b in str(fa.seq):
        B = b.upper()
        freq[B] += 1.0
    report(fa, freq)

