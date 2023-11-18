#!/bin/bash


#########################################################################################################
## SYNTAX:
##	trfBedUnitSummary.sh trf-file.tsv

## PREREQ:
## 	Convert to BED-like with: trfDatConverter.py -i trf-file.*dat > trf-file.tsv


## Created 2023-11-18.
##        Latest update on:       2023-11-18 
##        Other updates on:       2023-11-18 

#########################################################################################################


## INPUT
INBED=${1}

## HEADER
echo -e "#name\tunit_length\tunit_seq\tnum_arrays\ttotal_copy_number\ttotal_span\tavg_copy_per_array"

## BODY
awk '$1 !~/^#/ {a[$15]+=1 ; 		## number arrays
		b[$15]+=$5 ; 		## copy number (number units)
		c[$15]+=$3-$2 ;		## span
		d[$15]=$4 ;		## unit length ("period size"); "=" and not "+=" on purpose.
	}END{
	OFS="\t" ; 
	for (e in a) print 	d[e], 
				e, 
				a[e], 
				b[e], 
				c[e], 
				b[e]/a[e]}' ${INBED} | sort -k1,1n -k3,3nr -k4,4nr -k5,5nr -k6,6nr | awk 'OFS="\t" {print "unit_"NR, $1, $2,$3, $4, $5, $6}'

