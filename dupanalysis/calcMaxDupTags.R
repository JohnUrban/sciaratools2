## Copied this over from EmmuraldDir on 2021-07-12. Go there to make changes to master.

## USAGE: Rscript calcMaxDupTags.R 2*G N P
## EXAMPLE: Rscript calcMaxDupTags.R 1e6 1e5 0.001


## Read in arguments
args <- commandArgs(trailingOnly=TRUE)

## HELP
USAGE="\n	USAGE:	Rscript calcMaxDupTags.R 2*G N P  \n\n"
if(length(args)==0){cat(USAGE); quit()}


## PROCESS ARGS
genome_size <- as.numeric(args[1])

numTags <- as.numeric(args[2])

p <- as.numeric(args[3])


## FUNCTIONS
cal_max_dup_tags <- function(genome_size, numTags, p){
  for (x in 1:numTags){
    if (p > pbinom(x, numTags, 1/genome_size, lower.tail=FALSE)){
      return(x)
   }
  }
return(numTags)
}

## EXECUTE
maxTags <- cal_max_dup_tags(genome_size, numTags, p)

cmd <- paste0("echo ", maxTags)
system(cmd)
