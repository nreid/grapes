library(tidyverse)

# read in metadata table
meta <- read.table("../../../metadata/metadata.txt",stringsAsFactors=FALSE,sep="\t",comment.char="",header=TRUE)
rownames(meta) <- meta[,1]

# featurecounts fragments mapping within 2kb of restriction sites
co <- read.table("vvinifer_counts.txt",stringsAsFactors=FALSE,header=TRUE)
colnames(co) <- gsub(".bam","",colnames(co)) %>% gsub(".*\\.","",.)
# counts only
cot <- co[,-(1:6)]
cot <- as.matrix(cot)
# scaled by sum of mapped fragments
cots <- apply(cot,MAR=2,FUN=function(x){x/sum(x)})

# featurecounts summary
ss <- read.table("vvinifer_counts.txt.summary",stringsAsFactors=FALSE,header=TRUE,row.names=1)
colnames(ss) <- gsub(".bam","",colnames(ss)) %>% gsub(".*\\.","",.)

# samtools stats SN statistics
sn <- read.table("SN.txt",stringsAsFactors=FALSE,sep="\t",header=TRUE)

# make sure cot, ss, meta are in the same order as sn:
ss <- ss[,colnames(sn)]
cot <- cot[,colnames(sn)]
meta <- meta[colnames(sn),]

sitecounts <- rowSums(co[,-(1:6)])
sitecounts_l <- rowSums(co[,-(1:6)]) %>% log(.,10)
hist(sitecounts,breaks=500,xlim=c(0,100000))
hist(sitecounts_l,breaks=50)

# plot of raw sequences per sample
# 1M read cutoff
plot(sort(unlist(sn[1,])),ylab="raw sequences")
abline(h=1e6)

# use 1M sequences / individual as a cutoff. 
keepi <- unlist(sn[1,]) > 1e6

# pool 1 and 2 individuals
p1 <- grepl("Pool1",colnames(cot))
p2 <- grepl("Pool2",colnames(cot))

# histogram of total fragments per site
par(mfrow=c(1,3))

(rowSums(cot[,keepi])+1) %>% log(.,10) %>% hist(.,breaks=50,main="fragments per RAD site, log10")
abline(v=log(c(1000,80000),10))
(rowSums(cot[,keepi & p1])+1) %>% log(.,10) %>% hist(.,breaks=50,main="fragments per RAD site, log10")
abline(v=log(c(1000,40000),10))
(rowSums(cot[,keepi & p2])+1) %>% log(.,10) %>% hist(.,breaks=50,main="fragments per RAD site, log10")
abline(v=log(c(1000,40000),10))

# fragments per site for each pool in a scatter plot
	# many sites in pool1 not sequenced, or barely sequenced in pool2
	# these sites may contribute heavily to pool bias. 
	# keep sites within the box. 
plot(log(rowMeans(cot[,keepi & p2]+1),10),log(rowMeans(cot[,keepi & p1]+1),10),ylab="mean coverage per RAD site - pool 1",xlab="mean coverage per RAD site - pool 2")
abline(
	h=c(2,4),
	v=c(1.5,4))

keepl <- log(rowMeans(cot[,keepi & p2]+1),10) > 1.5 & 
		 log(rowMeans(cot[,keepi & p2]+1),10) < 4 &
		 log(rowMeans(cot[,keepi & p1]+1),10) > 2 &
		 log(rowMeans(cot[,keepi & p1]+1),10) < 4

# replot
plot(log(rowMeans(cot[keepl,keepi & p2]+1),10),log(rowMeans(cot[keepl,keepi & p1]+1),10),ylab="mean coverage per RAD site - pool 1",xlab="mean coverage per RAD site - pool 2")

# total fragments per site for retained individuals and loci
rowSums(cot[keepl,keepi]) %>% log(.,10) %>% hist(.,breaks=50)
rowSums(cot[keepl,keepi]) %>% hist(.,breaks=50)

# get rid of upper tail of coverage. 
	# after this step, 1363 rad sites retained
keepl <- keepl & rowSums(cot[,keepi]) < 75000


# histogram of missing data per site
	# even after filtering the worst individuals and sites
	# there is an average of 37% missing data
rowSums(cot[keepl,keepi]==0) %>% table() %>% plot()

# plot of missing data per individual, first sorted, then colored by pool
	# pool 1 individuals have uniformly lower missing data. 
	# pool 2 individuals have highly variable missing data. 
colSums(cot[keepl,keepi]==0) %>% sort() %>% plot()
colSums(cot[keepl,keepi]==0) %>% plot(.,col=factor(substring(names(which(keepi)),1,5)))

# x-y plot of total sequences vs missing data
	# not a strong correlation. 
plot(unlist(sn[1,keepi]),colSums(cot[keepl,keepi]==0),col=factor(substring(names(which(keepi)),1,5)))

# throw out individuals with more than 600 missing sites
	# only 16 individuals remain. 
keepi <- keepi & colSums(cot[keepl,]==0) < 600

# mds plot on scaled coverage
md <- (cots[keepl,keepi] + 1) %>% log() %>% t() %>% dist()
mds <- cmdscale(md)
plot(mds,col=factor(substring(rownames(mds),1,5)))

# variance in scaled coverage is way higher in pool 2
apply(cots[keepl,keepi],MAR=2,FUN=var) %>% plot()




