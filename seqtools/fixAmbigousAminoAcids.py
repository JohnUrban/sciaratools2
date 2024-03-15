#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from collections import defaultdict
from numpy.random import choice
parser = argparse.ArgumentParser(description="""

DESCRIPTION
        
    - Take a multi-fasta file (.fa file with >=1 entry) where sequences may contain non-canonical or ambiguous amino acid letters.
    - Output sanitized fasta that contains only canonical AAs.

    

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

    
    Recommended removal characters: * .
        - For example OorthoDB v10 had protein sequences that ended with a dot "."
        - OrthoDB v10 and v11 have "*" characters in the sequences.


    How to quickly check your FASTA for non-canonical and ambiguous letters...
    grep -i -E '[OUBZXJ]' ${FASTA} | grep -v ">" | less
    grep -v ">" ${FASTA} | grep -c -i -E '[OUBZXJ]'
    OR
    grep -i -E 'O|U|B|Z|X|J' ${FASTA} | grep -v ">" | less
    grep -v ">" ${FASTA} | grep -c -i -E 'O|U|B|Z|X|J'




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
parser.add_argument('--samplesize', '-n',
                   type=int, default=1000000, 
                   help='''Provide number of amino acids to sample to determine frequencies.
If 0, uniform probabilities are assumed.
If -1, the entire file will be used.
Default = 1 million. ''')
parser.add_argument('--pseudocount', '-c',
                   type=int, default=0, 
                   help='''You may want to provide a pseudocount (e.g. 1) if doing a small sample size.
Default is 0 to allow for correctly learning frequencies from sequences missing a certain base.''')
##parser.add_argument('--convertN', '-N', 
##                   action='store_true', default=False,
##                   help='''Convert Ns into ACGT. Default: False.''')
parser.add_argument('--upper', '-U', 
                   action='store_true', default=False,
                   help='''Return output sequence entirely upper-case instead of as is. Default: False.''')
parser.add_argument('--lower', '-L', 
                   action='store_true', default=False,
                   help='''Return output sequence entirely lower-case instead of as is. Default: False.''')
parser.add_argument('--keep_unexpected', '-K', 
                   type=str, default=False,
                   help='''When unanticipated letter is encountered, keep it as is instead of converting it using the X rule. Default: False.''')

parser.add_argument('--remove', '-r', 
                   type=str, default=False,
                   help='''Characters to remove when seen in protein sequence. For example, both * and . are seen in orthoDB protein sequences.
Recommended: --remove "*." (quotes needed).
Otherwise these characters will be treated as any other unanticipated character and replaced by the X rule.''')

parser.add_argument('--anticipate', '-A', 
                   type=str, default=False,
                   help='''To anticipate an otherwise un-anticipated character, follow this example for @:
-A @:ARND.
This will replace @ with A, R, N, or D.
This does not assume X and x are the same.
To anticipate both (or more than 1 letter), do:
-A @,!:ARND
That is comma-separate all that will be replace with same letters.
To replace another letter with a different rule, follow example:
-A @:ARND;!:KMFP
That is, keep things for same rule comma-separate and new rule starts after ;

Note that providing an already-antcipated letter (as defined way above) will result in that
letter-replacement rule being changed to what you define here.
Be careful.

Replacements need to be in the canonical amino acids -- all uppercase.

Note - if a replacement letter is not part of sample frequency, then
it will be given a count of 1 to avoid throwing an error.
This may be small compared to other counts.



Can use --pseudocount to avoid any ACGT being 0.
Can use --samplesize 0 so all are set to 1.
Note however that both will globally affect all replacement choices.

''')

args = parser.parse_args()





## functions
def probs(ltrs, freq):
    s = 0
    for b in ltrs:
        s += freq[b]
    p = []
    for b in ltrs:
        p.append( freq[b]/s )
    return p



safe_ltrs = "AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYy"
other_ltrs = "BbJjXxZzOoUu"


ctr = 0
pc = args.pseudocount
if args.samplesize == 0:
    freq = {"A":1.0,  "C":1.0,  "D":1.0,  "E":1.0,  "F":1.0,  "G":1.0,  "H":1.0,  "I":1.0,  "K":1.0,  "L":1.0,  "M":1.0,  "N":1.0,  "P":1.0,  "Q":1.0,  "R":1.0,  "S":1.0,  "T":1.0,  "V":1.0,  "W":1.0,  "Y":1.0}
else:
    if args.samplesize < 0:
        args.samplesize = float('inf')
    freq = {"A":0.0+pc,  "C":0.0+pc,  "D":0.0+pc,  "E":0.0+pc,  "F":0.0+pc,  "G":0.0+pc,  "H":0.0+pc,  "I":0.0+pc,  "K":0.0+pc,  "L":0.0+pc,  "M":0.0+pc,  "N":0.0+pc,  "P":0.0+pc,  "Q":0.0+pc,  "R":0.0+pc,  "S":0.0+pc,  "T":0.0+pc,  "V":0.0+pc,  "W":0.0+pc,  "Y":0.0+pc}
    for fa in SeqIO.parse(args.fasta, 'fasta'):
        for b in str(fa.seq):
            if ctr <= args.samplesize:
                ctr += 1
                B = b.upper()
                if B in "ACDEFGHIKLMNPQRSTVWY":
                    freq[B] += 1.0
            else:
                break


##    Anticipated ambiguous letters are in BbJjXxZz
##        - B converted to D or N
##        - J converted to I or L
##        - X converted to any canonical AA
##        - Z converted to E or Q
##
##    

translate = {}

translate['B'] = lambda : choice(['D','N'], p = probs('DN',freq))
translate['b'] = lambda : choice(['d','n'], p = probs('DN',freq))

translate['J'] = lambda : choice(['I','L'], p = probs('IL',freq))
translate['j'] = lambda : choice(['i','l'], p = probs('IL',freq))

translate['Z'] = lambda : choice(['E','Q'], p = probs('EQ',freq))
translate['z'] = lambda : choice(['e','q'], p = probs('EQ',freq))

translate['X'] = lambda : choice(["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"], p=probs('ACDEFGHIKLMNPQRSTVWY',freq))
translate['x'] = lambda : choice(["a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y"], p = probs('ACDEFGHIKLMNPQRSTVWY',freq))

translate['U'] = lambda : choice(["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"], p=probs('ACDEFGHIKLMNPQRSTVWY',freq))
translate['u'] = lambda : choice(["a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y"], p = probs('ACDEFGHIKLMNPQRSTVWY',freq))

translate['O'] = lambda : choice(["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"], p=probs('ACDEFGHIKLMNPQRSTVWY',freq))
translate['o'] = lambda : choice(["a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y"], p = probs('ACDEFGHIKLMNPQRSTVWY',freq))


if args.remove:
    for char in args.remove:
        translate[char] = lambda : ''


if args.anticipate:
    rulepairs = args.anticipate.split(";")
    for rule in rulepairs:
        keylist, replacements = rule.split(":")
        keys = keylist.split(",")
        replist = [e for e in replacements]
        for ltr in replist:
            try:
                assert freq[ltr] > 0
            except:
                freq[ltr] = 1.0
        for key in keys:
            translate[key] = lambda : choice(replist, p = probs(replacements,freq))

        
for fa in SeqIO.parse(args.fasta, 'fasta'):
    seq = str(fa.seq)
    newseq = ''
    for ltr in seq:
        if ltr not in safe_ltrs:
            try:
                newltr = translate[ltr]()
            except:
                #default to treating unantipated letters as X
                if not args.keep_unexpected:
                    newltr = translate['X']()
                else: #keep unexpexted
                    newltr = ltr
                
            newseq += newltr
        else:
            newseq += ltr
    if args.upper:
        newseq = newseq.upper()
    elif args.lower:
        newseq = newseq.lower()
    print(">"+fa.description)
    print(newseq)
        
        
