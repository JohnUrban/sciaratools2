#!/bin/bash


#########################################################################################################
## SYNTAX:
##	trfUnitSummaryStats.sh trf-file-unit-summary.tsv

## PREREQ:
## 	Convert to BED-like with: 	trfDatConverter.py -i trf-file.*dat > trf-file.tsv
##	Get Unit Summary:	  	trfBedUnitSummary.sh trf-file.tsv > trf-file-unit-summary.tsv

## Created 2023-11-18.
##        Latest update on:       2023-11-18
##        Other updates on:       2023-11-18

#########################################################################################################



## NOTE: header of input is : #name	unit_length	unit_seq	num_arrays	total_copy_number	total_span	avg_copy_per_array
## targeting 2=unit_length 4=num_arrays 5=total_copy_number 6=total_span 7=avgcopy per array


## INPUT
INSUMMARY=${1}


## SET UP
COLS=( 2 4 5 6 7 )
NAMES=( unit_length num_arrays total_copy_number total_span avg_unit_copy_per_array )




## OPEN PARENTHESIS FOR HEADER AND BODY
(
## HEADER
grep -v ^# ${INSUMMARY} | stats.py -k 2 -H -t -p 1,5,10,25,33,50,66,75,90,95,99 -A "." | head -n 1 | awk '{sub("string","METRIC"); gsub(",","\t"); print}'

## BODY
for i in {0..4} ; do
  grep -v ^# ${INSUMMARY} | stats.py -k ${COLS[$i]} -p 1,5,10,25,33,50,66,75,90,95,99 -t -A ${NAMES[$i]} | awk '{gsub(",","\t"); print}'
done 

## CLOSED PARENTHESIS FOR HEADER AND BODY :: PIPE ALL INTO AWK FOR REFORMATTING 
) | awk 'OFS="\t" {print $5,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$4,$6,$7,$8,$9,$10,$2,$3,$1}'



exit
