#!/usr/bin/env Rscript

library(stringr)
library(magrittr)

# This script finds restriction sites in a FASTA file
# It can be run as `cat genome.fasta | Rscript restriction.R >out.bed`
# It only considers 2 lines at at time, and prints results in BED format to the stdout

# for a non-palindromic cut site, do:
#	target <- c("AAAGGG","CCCTTT")

# here, restriction site sbfI
target <- c("CCTGCAGG")
tl <- nchar(target)


bl <- 0
basepos <- 0
baseseq <- c()
inseq <- FALSE

ln <- 1

f <- file("stdin")
open(f)
while(length(line <- readLines(f,n=1)) > 0) {

	# write(c(inseq, bl,basepos),stderr())

	if(grepl("^#",line)){inseq <- FALSE; ln <- ln + 1; next()}
	if(grepl("^$",line)){inseq <- FALSE; ln <- ln + 1; next()}
	# fasta header line
	if(grepl("^>",line)){
		sn <- substring(line,2)
		basepos <- 0
		baseseq <- c()
		inseq <- TRUE
		bl <- 0
		ln <- ln + 1; 
		next()
		}
	
	if(inseq){
		twolines <- paste(baseseq,line,sep="")

		cutsites <- str_locate_all(twolines,regex(target,ignore_case=TRUE)) %>% do.call(rbind,.)
		cutsites <- data.frame(cutsites)

		if(dim(cutsites)[1] == 0){
			basepos <- basepos + bl
			baseseq <- line
			bl <- nchar(baseseq)
			ln <- ln + 1
			next()
			}
			else{
				cutsites <- cutsites[cutsites[,1] > bl - tl[1] + 1,]
				if(dim(cutsites)[1] == 0){
					basepos <- basepos + bl
					baseseq <- line
					bl <- nchar(baseseq)
					ln <- ln + 1
					next()
					}
					else{

						cutsites <- cutsites + basepos
						# zero-index start position for bed format
						cutsites[,1] <- cutsites[,1] - 1
						cutsites <- cutsites[order(cutsites[,1]),]
						cutsites <- cbind(sn,cutsites)

						write.table(cutsites,stdout(),row.names=FALSE,col.names=FALSE,sep="\t")

						basepos <- basepos + bl
						baseseq <- line
						bl <- nchar(baseseq)
						ln <- ln + 1
						}
					}
	
	}

}
