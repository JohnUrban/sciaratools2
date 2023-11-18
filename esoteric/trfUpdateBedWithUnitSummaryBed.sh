#!/bin/bash


#########################################################################################################
## SYNTAX:
##	trfUpdateBedWithUnitSummaryBed.sh trf-file.tsv trf-file-unit-summary.tsv

## PREREQ:
## 	Convert to BED-like with: 	trfDatConverter.py -i trf-file.*dat > trf-file.tsv
##	Get Unit Summary:	  	trfBedUnitSummary.sh trf-file.tsv > trf-file-unit-summary.tsv

## Created 2023-11-18.
##        Latest update on:       2023-11-18
##        Other updates on:       2023-11-18

#########################################################################################################


## INPUT
INBED=${1}
INSUMMARY=${2}


## BODY
awk '$1 !~/^#/ {OFS="\t" ; 
		print 	$1,							## chr
			$2,							## start
			$3,							## end
			$15,							## unit_seq (will be translated to unit_name)
			$5,							## copy number
			".",							## strandless
			"unit_length:"$4";copy_number:"$5";unit_seq:"$15,	## PRIMARY INFO (will be BED name which can be viewed in IGV) ; Next line is remaining info.
			"Array_seq:"$16";Entropy:"$14";consensus_size:"$6";percent_matches:"$7";percent_indels:"$8";score:"$9";A:"$10";C:"$11";G:"$12";T:"$13	 ## OTHER INFORMATION
		}' ${INBED} | translateTable.py -i - -c 4 -d ${INSUMMARY} -k 3 -v 1 | awk 'OFS="\t" {a[$4]+=1 ; b[$1$4]+=1 ; c[$1]+=1 ; 
												print 	$1,							## Chr 
													$2, 							## Start
													$3, 							## End
													$4";unit_array:"a[$4]";"$1"_unit_array:"b[$1$4]";"$7, 	## Name with primary info
													$5, 							## Copy number (Score)
													$6, 							## Strandless placeholder
													"TR_array:"NR";"$1"_TR_array:"c[$1]";"$8		## Other Info
												}'
