#!/bin/bash

function help {
  echo "
  Usage:	md5EachSeqComparingTwoSeqs.sh fasta1 fasta2 seqnamesToCompare.txt

		Differs from md5EachSeqComparingTwoFiles.sh in following way:
		- the *Files.sh script extracts the fasta record and leaves the name and header as is, and sequence formatting can differ.
		- this *Seqs.sh script extracts the fasta record, formats the seq to be on a single line, makes both completely uppercase, and subtracts out the header.
			- thus this script will do better to ensure sequences are being compared, not records (that can differ with name/header, capitalizaton, and formatting).


  "
}

function errmsg {
  echo $@ 1>&2
}

if [ $# -eq 0 ]; then help ; exit ; fi

S=$( uname -s )
if [ $S == "Darwin" ]; then READ=greadlink ; else READ=readlink ; fi
R=${RANDOM}
WD=md5_wd_${R}
F1=$( ${READ} -f ${1} )
F2=$( ${READ} -f ${2} )
NAMES=$( ${READ} -f ${3} )

mkdir -p ${WD} && cd ${WD}

i=0
while read name ; do 
  let i++
  errmsg ${i} : ${name}
  extractFastxEntries.py -c $name -f ${F1} | fastaCaseMaker.py -f - | fastaFormatter.py -f - | grep -v ">" > f1.${name}.seq ; 
  extractFastxEntries.py -c $name -f ${F2} | fastaCaseMaker.py -f - | fastaFormatter.py -f - | grep -v ">" > f2.${name}.seq ; 
  A=( $( md5sum-lite f1.${name}.seq ) ) ; 
  B=( $( md5sum-lite f2.${name}.seq ) ) ;  
  echo -e "${A[1]}\t${B[1]}\t${A[0]}\t${B[0]}\t${name}" ; 
  rm f1.${name}.seq f2.${name}.seq
done <  ${NAMES}

## EXIT CLEANUP
cd ../ && rmdir ${WD}

