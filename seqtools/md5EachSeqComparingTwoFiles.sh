#!/bin/bash

function help {
  echo "
  Usage:	scriptname fasta1 fasta2 seqnamesToCompare.txt
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
  extractFastxEntries.py -c $name -f ${F1} > f1.${name}.fa ; 
  extractFastxEntries.py -c $name -f ${F2} > f2.${name}.fa ; 
  A=( $( md5sum-lite f1.${name}.fa ) ) ; 
  B=( $( md5sum-lite f2.${name}.fa ) ) ;  
  echo -e "${A[1]}\t${B[1]}\t${A[0]}\t${B[0]}\t${name}" ; 
  rm f1.${name}.fa f2.${name}.fa
done <  ${NAMES}

## EXIT CLEANUP
cd ../ && rmdir ${WD}

