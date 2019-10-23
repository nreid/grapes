library(tidyverse)
library(ape)
library(phytools)


# read in metadata table
meta <- read.table("../../../metadata/metadata.txt",header=TRUE,sep="\t",comment.char="",quote="",stringsAsFactors=FALSE)
rownames(meta) <- meta[,1]

# get vcf header line
# OS X line:
f <- pipe("gzcat ../results/freebayes/vvinifera_fb.vcf.gz | grep CHROM")
# LINUX line:
# f <- pipe("gzcat ../results/freebayes/vvinifera_fb.vcf.gz | grep CHROM")
h <- scan(f,what="character")

vcf <- read.table("../results/freebayes/vvinifera_fb.vcf.gz",stringsAsFactors=FALSE)
hcf <- read.table("../results/freebayes/vvinifera_fb_hap.vcf.gz",stringsAsFactors=FALSE)

colnames(vcf) <- h
colnames(hcf) <- h

keepi <- read.table("../../../metadata/retained_samples.txt",stringsAsFactors=FALSE) %>% unlist()

# keep high quality & biallelic variants
hq <- vcf[,6] > 50 & !grepl(",",vcf[,5])
vcf <- vcf[hq,]
hcf <- hcf[hq,]

# haploidized genotypes from filtered individuals only
hcfk <- hcf[,keepi] %>% as.matrix()
class(hcfk) <- "numeric"

# missing genotypes by site
rowSums(!is.na(hcfk)) %>% table() %>% plot()

# keep sites with more than X genotypes
keepl <- rowSums(!is.na(hcfk)) > 5 &
		 rowSums(hcfk,na.rm=TRUE) > 1 & 
		 rowSums(-1 * (hcfk - 1),na.rm=TRUE) > 1

vcf <- vcf[keepl,]
hcf <- hcf[keepl,]
hcfk <- hcfk[keepl,]

dd <- t(hcfk) %>% dist()
mds <- cmdscale(dd)
plot(mds,pch=20,col=factor(substring(rownames(mds),1,5)))

nj(dd) %>% midpoint.root() %>% plot()

