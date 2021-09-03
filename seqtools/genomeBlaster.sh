#!/usr/bin/env bash
set -e

function help {
  echo "
################################################################################
				GENOME BLASTER
################################################################################

	Usage: 
        
	genomeBlaster [s:q:S:Q:e:c:w:W:O:]
        	-s   Subject (target) genome.
		-q   Query genome.
		-S   Optional genome file for subject genome, used to speed up auto-word-size calculations.
		-Q   Optional genome file for query genome, used to speed up auto-word-size calculations.
	        -e   expected error rate. Provide decimal between 0-1. Default = 0.001. (1 error per 1000 bp, Default autowordsize 500).
		-c   Culling limit. Default = 1. Normally kept as is for this application.
		-w   Word size. Over-rides autowordsize calculation. Note that when using -c and/or -w, you might as well directly use blastimap2.2 or blastn, as this is mainly a convenience script.
		-W   Weight. Float between 0 and infinity, but typically <= 1. Used in auto wordsize calculation as W/error_rate. Default = 0.5. Not normally changed.
		-O   Other options to give BLAST. Provide in quotes.

	Step 1: if word_size not given, calculate smallest by.
		- if genome files not given, compute tmp files with fxSize.py (or faSize -detailed)
		- concatenate genome files, find shortest contig, divide that length by 2 and call it minlen
		- take expected error rate and compute 0.5*1/error_rate, calling this autowordsize
			- if error rate is 1%, then this is 0.5/0.01 = 50.
		- use word_size = min( minlen, autowordsize )
	Step 2: 
		- Use blastimap2.2 to align genomes and get PAF
			- larger word sizes will speed things up
			- so if you only care about the very longest scaffolds aligning, you can increase word_size to be pretty long
			- you can also just skip this wrapper script and use blastimap2.2 (or blastn) directly, specifying the long word_size and -culling_limit 1.

################################################################################
				GENOME BLASTER
################################################################################
"
}



##############################################################################
## DEFAULTS
##############################################################################
ERATE=0.001
WORDSIZE=-1
CULL=1
SEED=$RANDOM
WDIR=genomeblaster_tmp_wd_${SEED}
WEIGHT=0.5
MINLEN=-1
AUTO=-1

##############################################################################
## GET OPTS
##############################################################################
while getopts "s:q:S:Q:e:c:w:W:O:" arg; do
    case $arg in
        s) SFAS=${OPTARG};;
        q) QFAS=${OPTARG};;
        S) SGEN=${OPTARG};;
        Q) QGEN=${OPTARG};;
        e) ERATE=${OPTARG};;
        c) CULL=${OPTARG};;
        w) WORDSIZE=${OPTARG};;
        W) WEIGHT=${OPTARG};;
        O) OPTIONS="${OPTARG}";;
        *) help; exit;;
    esac
done


##############################################################################
## HELP CATCHALL
##############################################################################

if [ $# -eq 0 ]; then help; exit; fi		



## REPORT PARAMETERS - DEBUGGING
echo BEFORE AUTO CALCULATIONS 1>&2
for VAR in SFAS QFAS SGEN QGEN ERATE WEIGHT CULL MINLEN AUTO WORDSIZE ; do
  echo -e "${VAR}\t${!VAR}" 1>&2
done

##############################################################################
## STEP 1
##############################################################################

## REQUIRE SFAS AND QFAS BE PROVIDED
if [ -z $SFAS ] || [ -z $QFAS ] ; then echo "Subject and query fasta files are required arguments." ; help; exit; fi

## IF NEED TO CALCULATE WORD SIZE; THEN DO STUFF
if [ ${WORDSIZE} == -1 ]; then 
	## IF NO GENOME FILE FOR SUBJECT, THEN GENERATE
	if [ -z $SGEN ]; then mkdir -p ${WDIR} ; faSize -detailed ${SFAS} > ${WDIR}/subject.genome ; SGEN=${WDIR}/subject.genome; fi

	## IF NO GENOME FILE FOR SUBJECT, THEN GENERATE
	if [ -z $QGEN ]; then mkdir -p ${WDIR} ; faSize -detailed ${QFAS} > ${WDIR}/query.genome ;  QGEN=${WDIR}/query.genome; fi
	
	## GET MIN CONTIG LEN AND DIVIDE BY 2
	MINLEN=$( cat $SGEN $QGEN | sort -k2,2n | head -n 1 | awk '{print $2/2}' )

	## GET AUTOWORDSIZE
	AUTO=$( echo ${WEIGHT}/${ERATE} | bc -l )

	## GET MINIMUM WORDSIZE
	WORDSIZE=$( echo -e "${MINLEN}\n${AUTO}" | sort -n | head -n 1 | awk '{print int($1)}' )

	## CLEANUP
	if [ -d ${WDIR} ]; then echo $WDIR 1>&2 ; ls $WDIR 1>&2 ; rm -r ${WDIR} ; fi
fi	

## REPORT PARAMETERS
echo AFTER AUTO CALCULATIONS 1>&2
for VAR in SFAS QFAS SGEN QGEN ERATE WEIGHT CULL MINLEN AUTO WORDSIZE ; do
  echo -e "${VAR}\t${!VAR}" 1>&2
done


## RUN BLASTIMAP
echo "blastimap2.2 ${SFAS} ${QFAS} -word_size ${WORDSIZE} -culling_limit ${CULL} ${OPTIONS}"
blastimap2.2 ${SFAS} ${QFAS} -word_size ${WORDSIZE} -culling_limit ${CULL} ${OPTIONS}



## NOTES TO SELF
## 1000 bp word size for alns -- means contigs shorter than 1000 wont be aligned... I think.
# blastimap2.2 ../data/canu_primary-chromosome-correct-oritentation.fasta ../data/canu_all-chromosome-correct-oritentation.fasta -word_size 1000 -culling_limit 1 > del.paf

## A formula for making a word size depending on how much agreement you expect
#- If you expect 99% similarity, 1% divergence
#	- Then the longest you can go is 99 bp, but only 1 99 bp stretch per 100 bp will have no error, so you'd want more seeds...
#		- the question is how many seeds do you want per 100 bp in this scenario? 99 bp would be 1/100 or 10/1000, or 100/10e3 and so on.. is that enough seeds?
#		- can simply do X/2
#		- so for 99% similarity expectation, 99/2 = 49 or 50
#		- but 99.999 would be the same...
#		- so really you care about 1/0.99 and 1/0.99999	= 1.01 and 1.00001
#		- no, you care about 1/0.01 and 1/0.00001 = 100 and 10,000
#		- so the forumula could be:
#			- 0.5*1/error_rate = suggested word_size
