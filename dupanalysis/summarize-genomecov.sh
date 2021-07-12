#!/bin/bash

function help {
echo "
	Usage: 
		summarize-genomecov.sh SHARED_PREFIX [K]
			SHARED_PREFIX = Shared prefix for bedtools genomecov files in same dir. See below.
			K = Max number of 5 prime ends or max coverage over each base to "keep" at each position.

	Upstream requirements:
		1. For a given BAM file, do:
			BASE=$( basename $f .bam );  
			bedtools genomecov -5 -ibam ${BAM} -g ${G} > ${OUTDIR}/${BASE}.txt & 
			bedtools genomecov -strand "+" -5 -ibam ${f} -g ${G} > ${OUTDIR}/${BASE}.pos.txt & 
	  		bedtools genomecov -strand "-" -5 -ibam ${f} -g ${G} > ${OUTDIR}/${BASE}.neg.txt & 

		   The script as currently written actually only uses the pos and neg files above.

		2. Determine K values to analyze.
			Use: calcMaxDupTags.R 
				Usage: 
					Rscript calcMaxDupTags.R 2*G N P
						INPUTS
					 	G = Genome length (multiply by 2 to account for both strands)
						N = number of reads or pairs
						P = p-value to define cutoff between unsurprising and surprising.
										
						OUTPUT = K:
						K = the max number of reads at a site considered as unsurprising.

	    		Smaller p-value cutoffs (closer to 0) gives larger K.
			Bigger p-value cutoffs (closer to 1) gives smaller K.
			This has confusing effect where p-value cutoffs closer to 1 (bigger) is more conservative.
				i.e. Considers fewer reads "real", more as "duplicates".
			So avoid urge to use a really small p-value cutoff as a way to be conservative.
			Recommended p-value cutoff range:
			0.1 to 0.00001

			Using the same 2*G and P values, try altering N value to see K values for enriched (or depleted) regions.
			For example, you may have an average of 10-fold enriched peak regions.
				Use 10*N, to get an idea of K for those regions.
			Or you may have 2-fold depleted regions:
				Use N/2 to get an idea for K in those regions.
												    

	Summary:
		A common step during genomic analysis, after mapping reads, is to remove redundant reads.
		This is an attempt to remove PCR (and optical) artifacts.
		Various tools determine whether reads should be labeled as duplicate or not (and possibly removed).
		A shortcoming of most tools is applying 1 rule to every position.
		For example, Picard labels all except 1 putative duplicates (per position) as duplicate (and optionally removes).
		As another example, the common peak calling pipeline (MACS2) offers its own redundant read removal.
			MACS2 understands that more duplicate reads are expected for larger libraries (given the same genome size).
			It offers the options of removing all (keeping 0), keeping 1 (like Picard), or keeping K.
			K is some number determined to be unsurprising given the number of reads and genome size.
			So K can be thought of as arising from a genomic average.
			This single rule fails in the following ways:
				Depleted regions should have a smaller K than the genomic average.
				Enriched regions should have a larger K than the genomic average.
		This tool summarizes the output of BEDtools genomecov of 5 prime ends.
		It gives a different view of the duplication rates.
		The emphasis is on how many sites in the genome are affected by different putative duplications.
		It also calculates duplication rates (per base and per read) given you'd be fine leaving K at each site.
			So, this therefore still only uses a single rule for all sites.
			However, try different values of K to get a better idea.
				The first value of K to try would be the MACS2 formula for the number of reads you have and genome size.
				Another K value to try is using the MACS2 formula with your number of reads multiplied by an average or max fold-enrichment of peaks / enriched regions.
					This shows how many reads would be acceptable to leave in an enriched region.

	Full pipeline:
		1. bowtie2 to map and get BAM
		2. BEDtools to get coverage of 5prime ends
		3. calcMaxDupTags for K values to analyze
		4. summarize-genomecov given BEDtools files and K values.
		5. Iterate as need be.

" ;
}

if [ -$# -eq 0 ]; then help ; exit; fi


F=$1
K=5
if [ $# -gt 1 ]; then K=$2; fi

B=$( basename $F .txt )
D=$( dirname $F )
PRE=${D}/${B}

for C in ${PRE}.pos.txt ${PRE}.neg.txt; do
  grep genome ${C} | head -n 2 ; 
  grep genome ${C} | awk '$2>1 {s+=$5}END{print "genome\t2+\t-\t-\t"s}'
  grep genome ${C} | awk '$2>1 {s+=$2*$3}END{print "genome\t2+\t-\tNreads\t"s}'
  grep genome ${C} | awk '$2>1 {s+=($2-1)*($3)}END{print "genome\t2+\t-\tNreadsLeave1\t"s}'
  echo
  grep genome ${C} | awk -v "K=$K" '$2>K {s+=$5}END{print "genome\t"K+1"+\t-\t-\t"s}'
  grep genome ${C} | awk -v "K=$K" '$2>K {s+=$2*$3}END{print "genome\t"K+1"+\t-\tNreads\t"s}'
  grep genome ${C} | awk -v "K=$K" '$2>K {s+=($2-K)*($3)}END{print "genome\t"K+1"+\t-\tNreadsLeaveK\t"s}'
  echo
  echo
done
