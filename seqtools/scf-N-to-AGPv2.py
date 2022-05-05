#!/usr/bin/env python3

import argparse, sys
import re
from Bio import SeqIO

# Parse command-line arguments
parser = argparse.ArgumentParser(description='''
This is version 2.
- This update aimed at writing AGP files with gaps of both known and unknown length.
- The previous version only wrote AGPs with known length.
- I am leaving that version as is before the update b/c the update may require bulk changes to code.
- This will approach this in two ways:
    - 1. Use a specified length to identify gaps of unknown length -- this is 100 bp by NCBI standards.
    - 2. Interpret known and unknown gap characters in scaffolds (such as N and n).
- There are no use cases with "recognition" sequences in gaps of unknown size.
    - These are usually in optical map scaffolds that create gaps of known sizes.
    - Gaps of unknown size (U gaps) are completely unknown, often from HiC these days.
    - This script will, for now, have no option to search for rec sites inside U gaps.
    - This script will keep tab of Ngap and Ugap statistics.
    - It will raise warning to user if Ugaps not all one size, as that is frowned upon by NCBI.
- The script is automatically detecting both character types for Ngaps and Ugaps, but either can be turned off.

===============================================================================================
Description:

INPUT:
    - Takes in Fasta of scaffolds with N gaps between contigs with the goal of creating NCBI-acceptable AGP/ctg/scf files for it.

OUTPUT(S):
    - AGP
    - contigs fasta


DEPRECATED OUTPUTS:
    - Scaffolds fasta file - as it would be reproduced using the AGP and Contigs output files:
        - Use contigs-and-AGP-to-scaffolds.py to re-build scaffolds, or build them differently according to AGP instructions output here.
        - can use --outputscaffolds to get input scaffs returned... perhaps just as a sanity check.
        
    - NOTE on Contigs fasta file.
        - Can use AGP file to construct BED file to use with fastaFromBed or extractFastxEntries.py --bed, to use on input scaffolds (unless instructed changes to scaff construction in AGP).

OUTPUTS TODO 
    - TODO (OR UPDATE scf-N-to-BED.py)
        - also output 5-col gap bed file: chr, start, end, gaptype, gaptype_idx, gaplength
        - gaptype is "N" or "U"
        - gaptype_idx corresponds to the number of that gaptype seen including this gap (in order)
        - gaplength will allow quick stats with stats.py
        - actually this doesn't have to be limited to the gaps, it can contain the contig coords as well
        - it could basically be a BED-file version of the AGP file that maps/annotates/partitions the scaffolds.

For now, it is hard-coded to assume you're breaking apart a BioNano Map scaffolded assembly.

This script is similar to scf-N-to-BED.py (copied from there for starters).

Use '-' or 'stdin' if coming from stdin.
Use -l/--length to control minimum gap length (default 25).
Use -r/--recognition to allow regex to expand gaps through recognition sequence islands dispersed in N-gaps (e.g. from optical maps).
When using -r provide only the recognition sequence. For now it is case-sensitive.
By default N-gaps are simply where Ns start to the first non-N.

Example recognition sequences:
BssSI: CACGAG
BspQI: GCTCTTC

AGP coordinates are 1-based inclusive.

AGP entry for contig:
name start end part_number component_type name_id contig_start contig_end orientation
component_type = W
orientation = +


AGP entry for gap:
name start end part_number component_type gapLen gapType linkage evidence
component_type = N
gapType = scaffold
linkage = yes
evidence = map


Example of scaffold with 2 contigs separated by one gap from map-evidence:
trial_scaffold_1        1       9       1       W       trial_scaffold_1_1      1       9       +
trial_scaffold_1        10      33      2       N       24      scaffold        yes     map
trial_scaffold_1        34      39      3       W       trial_scaffold_1_3      1       6       +







=======================================================================
Historical context for version 1:
I started using a similar tool (fasta2agp.pl) found here:
    wget --no-check-certificate --content-disposition https://gist.githubusercontent.com/IdoBar/8a324fff97c6a866ea7e0cbcec41d570/raw/7faa9e342d62355b11c79f114ec6367e1de9525e/fasta2agp.pl


- I found that fasta2agp.pl can only filter by scaffold size, but
        - can't discriminate gap sizes (breaks at 1 bp)
        - can't mow over recseqs (so treats as multiple gaps)
            - BioNano scaffolding has huge gaps that can have multiple short (6 bp) recognition sequences within them
            - This can be useful in some contexts
            - It may not be useful to treat them all as contigs though...
            
- I created my own tool to do the same thing as fasta2agp.pl with the additional features I needed.
- It also default outputs what I needed.
        - It is called:
            - scf-N-to-AGP.py (this script)
        - It follows fasta2agp.pl behavior on contigs that do not start with N
                - If contig starts with N
                        - fasta2agp.pl ignores it as if it were not there and all coordinates are shifted
                        - scf-N-to-AGP.py gives a warning, but includes it as a "part" of the object
        - scf-N-to-AGP.py is able to mow over the recseqs and only break contigs at gap sizes exceeding some minlength.

        

''', formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("fasta")
parser.add_argument("-l", "--mingaplength", type=int, default=25, help='''Do not report as gap if less than this length. Default=25.''')
parser.add_argument("-r", "--recognition", type=str, default=False, help='''Allow N-gaps to encompass occurences of given recognition sequence. Default=False.''')
parser.add_argument("-G", "--ngapchar", type=str, default='N', help='''The single letter gap character for gaps of known length (Ngaps). This will produce N components in AGP, unless --ugaplen is used, in which case gaps found to be that specified length will be assumed to be U components. Default =  N.''')
parser.add_argument("-g", "--ugapchar", type=str, default='n', help='''The single letter gap character for gaps of unknown length (Ugaps; will produce U components in AGP). Default =  n.''')
parser.add_argument("-u", "--ugaplen", type=int, default=False, help='''Not used by default. If set, gaps found that equal the length given will be assumed to gaps of unknown length (U components in AGP).''')
#parser.add_argument("-c", "--ugaplenout", type=int, default=False, help='''If .''')
parser.add_argument("-P", "--prefix", type=str, default='', help='''Add this prefix in front of all sequence names.''')
parser.add_argument("-o", "--outprefix", type=str, default='', help='''Add this prefix in front of all output filenames. Can be path to directory or path + prefix as well.''')
parser.add_argument("-N", "--Ngapsonly", action='store_true', default=False, help='''Tells script to only anticipate Ngap characters.. Writes N components only.''')
parser.add_argument("-U", "--Ugapsonly", action='store_true', default=False, help='''Tells script to only anticipate Ugap characters.. Writes U components only.''')
caseparser = parser.add_mutually_exclusive_group()
caseparser.add_argument("-B", "--toupper", action='store_true', default=False, help='''Tells script to convert all characters in sequences to uppercase before searching.The use case is if you have only gaps of known length, but suspect there may be mixtures of upper- (N) and lower-case (n) gapchars from upstream processing (e.g. masking, etc). The script will see all as upper-case gapchar (N).''')
caseparser.add_argument("-S", "--tolower", action='store_true', default=False, help='''Tells script to convert all characters in sequences to lowercase before searching. The use case is if you have only gaps of unknown length, but suspect there may be mixtures of upper- (N) and lower-case (n) gapchars from upstream processing (e.g. masking, etc). The script will see all as lower-case gapchar (n).''')
parser.add_argument("-Z", "--tigSeqOptOut", type=int, default=1, help='''Tells script how to write contig sequences. There are two options (use the integer corresponding to desired option). 1 = the original input sequence before any processing. 2 = the sequence after processing (e.g. --toupper or --tolower). This option does not affect gap characters. Default = 1.''')
#parser.add_argument("-X", "--ngapcharout", type=str, default='N', help='''When re-writing the scaffolds, this tells the script what character to use for gaps of known length (Ngaps). Default = N.''')
#parser.add_argument("-x", "--ugapcharout", type=str, default='n', help='''When re-writing the scaffolds, this tells the script what character to use for gaps of unknown length (Ugaps). Default = n. A good option could be "X" or "x" to avoid upper/lower case issues.''')
parser.add_argument("--outputscaffolds", action='store_true', default=False, help='''Will return the input scaffolds in a new file. Somewhat pointless, except for a sanity check. In the future, this may output scaffolds updated according to instructions, but for now it is recommended to use the output contigs with the AGP file in the contigs-and-AGP-to-scaffolds.py script.''')

args = parser.parse_args()



def get_contig_info(record, gap, scfstart, part, prefix=''):
    ## TODO: where are scfstart and scfend come from.... is everything ok here?!
    
    scfend = gap[1].start() # No need to subtract 1 as this is 0-based (need to add 1 for real gap start)
    ctgend = scfend - scfstart + 1
    return (str(e) for e in (prefix+str(record.id),
                             scfstart,
                             scfend,
                             part,
                             'W',
                             prefix+str(record.id)+'_'+str(part),
                             1,
                             ctgend,
                             '+'))


def get_last_contig_info(record, gap, scflength, scfstart, part, prefix=''):
    ## ... "gap" variable does not seem to be used here. TODO - remove...
    scfend = scflength
    ctgend = scfend - scfstart + 1
    return (str(e) for e in (prefix+str(record.id),
                             scfstart,
                             scfend,
                             part,
                             'W',
                             prefix+str(record.id)+'_'+str(part),
                             1,
                             ctgend,
                             '+'))


def get_gap_info(record, gap, part, gapType='scaffold', linkage='yes', prefix=''):
    ''' For now Ngaps assumed to be map and Ugaps assumed to be hi-c paired-ends.'''
    gapLen = getgaplen(gap[1])
    gdict = {"N":"map",
             "U":"paired-ends"}
    if gap[0] not in "NU":
        # yell and quit
        error_message("Unanticipated gap type encountered. I only accept N and U.")
    ## TODO: allow ugaplenout to be specified, but ALSO NEED TO CHANGE IN SCAFFOLDS OUTPUT... unless I just scrap the entire idea of outputting tigs and scaffs... or scaffs..
    return (str(e) for e in ( prefix+str(record.id),
                                  gap[1].start()+1,
                                  gap[1].end(),
                                  part,
                                  gap[0],
                                  gapLen,
                                  gapType,
                                  linkage,
                                  gdict[ gap[0] ]
                                  ))


##def get_Ngap_info(record, gap, part, gapType='scaffold', linkage='yes', evidence='map', prefix=''):
##    return get_gap_info(record=record, gap=gap, part=part, gapType=gapType, linkage=linkage, evidence=evidence, componentType="N", prefix=prefix)
##
##def get_Ugap_info(record, gap, part, gapType='scaffold', linkage='yes', evidence='paired-ends', prefix=''):
##    return get_gap_info(record=record, gap=gap, part=part, gapType=gapType, linkage=linkage, evidence=evidence, componentType="U", prefix=prefix)

def get_file_name(x,pre, sfx):
    if not pre:
        pre = x
    else:
        pre = pre + '.' + x
    name = pre + '.' + sfx
    sys.stderr.write('Will write to ' + name + ' ....\n')
    return name

def error_message(msg=""):
    sys.stderr.write("ERROR: " + msg + "\n")
    sys.stderr.write("You can't fire me if I QUIT.\n")
    quit()
    
def returnseq(seq, toupper=False, tolower=False, returnOpt=1):
    '''
    seq = string
    toupper = optionally make all letters uppercase ONLY IF tigSeqOptOut=2
    tolower = optionally make all letters lowercase ONLY IF tigSeqOptOut=2
    tigSeqOptOut = whether to output sequence string as it was given (1, as-is) or processed according to case instructions if any given, else as-is (2, possibly-processed).
    '''
    assert not (toupper and tolower)
    if returnOpt == 1: ## return seq as is.
        return seq
    elif returnOpt == 2:
        if toupper:
            return seq.upper()
        elif tolower:
            return seq.lower()
        else: ## no processing was specified; return seq as is.
            return seq 
    else:
        # yell and quit
        error_message("Encountered ineligible tigSeqOptOut option. Only accepts 1 or 2.")
        
def append_sequence_to_fasta(fh, record, prefix='', sep='', part='', start=-1, end=-1, toupper=False, tolower=False, tigSeqOptOut=1):
    '''
    fh = file connect in write mode
    record = Bio fasta record
    prefix = optional string to add in front of record name
    sep = separator/delimiter between record.id and part
    start = optional start position in record (-1 tells to use whole record)
    end = optional end position in record (-1 tells to use whole record)
    toupper = optionally make all letters uppercase ONLY IF tigSeqOptOut=2
    tolower = optionally make all letters lowercase ONLY IF tigSeqOptOut=2
    tigSeqOptOut = whether to output sequence as is (1) or processed according to case instructions (2).

    tigSeqOptOut obviously seems redundant in the context of this function, but is needed in the larger context of this entire pyscript.
    '''
    # Name
    name = '>'+str(prefix)+str(record.id)+str(sep)+str(part)+'\n'
    
    # Sequence
    baseseq = returnseq(str(record.seq), toupper, tolower, tigSeqOptOut)
    if start == -1 or end == -1:
        # Return whole record (probably to scaffolds fasta)
        sequence = returnseq(str(record.seq), toupper, tolower, tigSeqOptOut) + '\n'  #returnOpt is tigSeqOptOut here
    else:
        # Return slice of record (probably to scaffolds fasta) (. . . did I mean probably to contigs fasta here?)
        sequence = returnseq(str(record.seq)[start:end], toupper, tolower, tigSeqOptOut) + '\n' #returnOpt is tigSeqOptOut here
        
    ## Write
    fh.write(name)
    fh.write(sequence)


def getgaplen(gap):
    return gap.end()-gap.start()


    
# OPEN CONNECTION TO FILE/STDIN (( INPUT SCAFFOLDS IN FASTA FORMAT ))
handle = sys.stdin if args.fasta in ('-','stdin') else open(args.fasta)



# DEFINE REGEX for NGAPS 
#N = args.gapchar + '+' ## 'N+'
Npattern = args.ngapchar + '+' ## 'N+'
if args.recognition:
    Npattern = N + '((' + args.recognition + '){0,1}' + N + ')*' ## 1+ Ns and optionally 0 or more (recseq, Ns) 


# DEFINE REGEX for UGAPS 
Upattern = args.ugapchar + '+' ## 'n+'

# Shorten some variable names
P = args.prefix
S = args.tolower
B = args.toupper
T = args.tigSeqOptOut

# Open contig, fasta, and AGP files for writing (W MODE)
if args.outputscaffolds:
    scaffolds_fh = open( get_file_name( 'scaffolds', args.outprefix, 'fasta'), 'w')
contigs_fh = open( get_file_name( 'contigs', args.outprefix, 'fasta'), 'w')
AGP_fh = open( get_file_name( 'AGP', args.outprefix, 'agp'), 'w')



# LOOP OVER FASTA
for record in SeqIO.parse(handle, "fasta"):
    
    if args.outputscaffolds:
        # Write out scaffold
        # TODO - 5/30/21 - the scaffold sequence needs to be re-built given the instructions of ugaplenout, ngapcharout, ugapcharout, tigSeqOptOut, etc...
        #
        append_sequence_to_fasta(scaffolds_fh,
                                 record,
                                 prefix=P,
                                 toupper=args.toupper,
                                 tolower=args.tolower,
                                 tigSeqOptOut=args.tigSeqOptOut
                                 ) ## tigSeqOptOut must be set to args.tigSeqOptOut here
    
    # Length of scaffold
    scflength = len(record.seq)

    # The search-sequence may optionally be processed -- if there is no processing to be done, the input string is returned as-is
    sequence2search = returnseq(str(record.seq),
                                args.toupper,
                                args.tolower,
                                returnOpt=2
                                ) ## returnOpt must be set to option 2 here to follow instructions for args.toupper/tolower if given.

    # Putative NGap locations -- Will search for Ngapchar, but if unknown gaps are a specified length, any found at that length will be re-classified as Ugaps
    if not args.Ugapsonly:
        allNgaps = [['N',e] for e in re.finditer(Npattern, sequence2search)]
    else:
        allNgaps = []

    # Ugap locations as determined by ugapchar. Ngaps of a specified length may be re-classified as U later. So, this is not the final number.
    if not args.Ngapsonly:
        allUgaps = [['U',e] for e in re.finditer(Upattern, sequence2search)]
    else:
        allUgaps = []

    # Possible re-classification of Ngaps as Ugaps given a specified length
    nReclassified = 0
    if args.ugaplen is not False and not args.Ugapsonly:
        for i in range(len(allNgaps)):
            gap = allNgaps[i][1]
            if getgaplen(gap) == args.ugaplen:  ## args.ugaplen is "False" by default.
                nReclassified += 1
                allNgaps[i][0] = 'U'

    # Combine and sort gaps
    allgaps = allNgaps + allUgaps
    allgaps.sort(key = lambda x: x[1].start())
    
    # Number of gaps found
    numUgaps = len(allUgaps) + nReclassified
    numNgaps = len(allNgaps) - nReclassified
    numgaps = len(allgaps)

    # Sanity checks
    assert numgaps == numUgaps + numNgaps
    assert numgaps == len(allUgaps) + len(allNgaps)
    
    # Process scaffold accordingly
    if numgaps == 0:
        # Single contig - write AGP 
        outline = ("\t").join([str(e) for e in (P+str(record.id), 1, scflength, 1, 'W', P+str(record.id)+'_1', 1, scflength, '+')])
        AGP_fh.write( outline + '\n')
        # Write contig 
        append_sequence_to_fasta(contigs_fh,
                                 record,
                                 prefix=P,
                                 sep='_',
                                 part=1,
                                 toupper=args.toupper,
                                 tolower=args.tolower,
                                 tigSeqOptOut=args.tigSeqOptOut)

    elif numgaps > 0: # 1 or more gaps
        
        # Describe initial contig
        scfstart = 1
        part = 1
        ended_on_gap = False
        
        for gap in allgaps:
            # If current gap is not minimum gap length, do not report it, do not update current scfstart variable, just continue to next gap
            if getgaplen(gap[1])  < args.mingaplength:
                continue
            
            # Warn about leading Ns if present, report gap, and move on
            if gap[1].start() == 0:
                sys.stderr.write("Warning: "+str(record.id)+" started with N. The leading N-gap will be reported as the first part of this scaffold. NCBI does not accept this though.\n")

                # Report gap info
                ## TODO: 5/30/21 - gap processing needs to be updated to deal with gap-list-elements and N vs U
                AGP_fh.write( ("\t").join(
                    get_gap_info( record,
                                  gap,
                                  part,
                                  prefix=P)) + '\n' )

                # Update scfstart s.t. it starts after this gap
                ## TODO: 5/30/21 - this handling may have to change according to ugaplenout....
                scfstart = gap[1].end()+1

                #Update part number for next item
                part += 1
                continue

            

            # Process previous contig (and update part number):
            #   Since iterating over gaps, contigs precede current gap starting at most recent scfstart variable and ending at gap

            # Write AGP
            ## TODO: 5/30/21 - gap processing needs to be updated to deal with gap-list-elements and N vs U
            ## TODO: 5/30/21 - this handling may have to change according to ugaplenout....
            AGP_fh.write( ("\t").join(
                get_contig_info( record,
                                 gap,
                                 scfstart,
                                 part,
                                 prefix=P) ) + '\n' )

            # Write contig
            ## TODO: 5/30/21 - this handling may have to change according to ugaplenout....
            append_sequence_to_fasta(contigs_fh,
                                     record,
                                     prefix=P,
                                     sep='_',
                                     part=part,
                                     start=scfstart-1,
                                     end=gap[1].start(),
                                     toupper=args.toupper,
                                     tolower=args.tolower,
                                     tigSeqOptOut=args.tigSeqOptOut)
            # Update part
            part += 1


            # Process current gap, update scfstart for next contig and part number
            ## TODO: 5/30/21 - gap processing needs to be updated to deal with gap-list-elements and N vs U
            ## TODO: 5/30/21 - this handling may have to change according to ugaplenout....
            AGP_fh.write( ("\t").join(
                get_gap_info(record,
                             gap,
                             part,
                             prefix=P)) + '\n' )
            
            scfstart = gap[1].end()+1
            part += 1
            
            # Warn about trailing Ns if present, report gap (above), and break
            if gap[1].end() == scflength:
                ended_on_gap = True
                sys.stderr.write("Warning: "+str(record.id)+" ended with N. The trailing N-gap will be reported as the last part of this scaffold. NCBI does not accept this though.\n")
                break



        # Process last contig (if it did not end on a gap)
        if not ended_on_gap:
            ## TODO: 5/30/21 - gap processing needs to be updated to deal with gap-list-elements and N vs U
            ## TODO: 5/30/21 - this handling may have to change according to ugaplenout....
            # Write AGP
            AGP_fh.write( ("\t").join(
                get_last_contig_info(record,
                                     gap,
                                     scflength,
                                     scfstart,
                                     part,
                                     prefix=P) ) + '\n' )
            # Write contig
            append_sequence_to_fasta(contigs_fh,
                                     record,
                                     prefix=P,
                                     sep='_',
                                     part=part,
                                     start=scfstart-1,
                                     end=scflength,
                                     toupper=args.toupper,
                                     tolower=args.tolower,
                                     tigSeqOptOut=args.tigSeqOptOut)
                


# Close input file            
if args.fasta not in ('-','stdin'):
    handle.close()


# Close output files
if args.outputscaffolds:
    scaffolds_fh.close()
contigs_fh.close()
AGP_fh.close()



# END


