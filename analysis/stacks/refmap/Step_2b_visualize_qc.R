library(tidyverse)

# this script is used to explore coverage among
# rad sites and among individuals. 

# in the end, a bed file is output giving a set of RAD sites
# to call variants on. 

# assumes you're running the script interactively from the directory you found it in. 

# read in metadata table
meta <- read.table("../../../metadata/metadata.txt",stringsAsFactors=FALSE,sep="\t",comment.char="",header=TRUE)
rownames(meta) <- meta[,1]

# featurecounts fragments mapping within 2kb of restriction sites
co <- read.table("../results/align_stats/vvinifer_counts.txt",stringsAsFactors=FALSE,header=TRUE)
colnames(co) <- gsub(".bam","",colnames(co)) %>% gsub(".*\\.","",.)
# matrix with counts only
cot <- co[,-(1:6)]
cot <- as.matrix(cot)
# scaled by sum of mapped fragments
cots <- apply(cot,MAR=2,FUN=function(x){x/sum(x)})

# featurecounts summary
ss <- read.table("../results/align_stats/vvinifer_counts.txt.summary",stringsAsFactors=FALSE,header=TRUE,row.names=1)
colnames(ss) <- gsub(".bam","",colnames(ss)) %>% gsub(".*\\.","",.)

# samtools stats SN statistics
sn <- read.table("../results/align_stats/SN.txt",stringsAsFactors=FALSE,sep="\t",header=TRUE)

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

# the coverage distribution is so screwed up
# that I'll use # of missing rad sites < 1600 as a cutoff
# instead of total coverage. 
# 20 individuals retained
keepi <- colSums(cot==0) < 1600

# pool 1 and 2 individuals
p1 <- grepl("Pool1",colnames(cot))
p2 <- grepl("Pool2",colnames(cot))


# pooling is highly unequal. 
	# to examine this, calculate the
	# effective number of samples
	# per https://en.wikipedia.org/wiki/Diversity_index
	# maybe not perfectly relevant to downstream analyses, 
	# but a good measure of inequality of pooling. 


# p is a vector of proportions
# q is the hill number order
qD <- function(p,q){

	if(q == 1){
		out <- exp(-sum(p * log(p)))
	}
	if(q > 1){
		out <- (sum(p^q))^(1/(1-q))
	}

	return(out)
}

# effectively, there are only 11.5 samples in pool 1
qD(unlist(sn[1,p1]/sum(sn[1,p1])),1)

# effectively, there are only 18 samples in pool 2
qD(unlist(sn[1,p2]/sum(sn[1,p2])),1)



# histogram of total fragments per site
par(mfrow=c(1,3))

(rowSums(cot[,keepi])+1) %>% log(.,10) %>% hist(.,breaks=50,main="fragments per RAD site, log10")
(rowSums(cot[,keepi & p1])+1) %>% log(.,10) %>% hist(.,breaks=50,main="fragments per RAD site, log10")
(rowSums(cot[,keepi & p2])+1) %>% log(.,10) %>% hist(.,breaks=50,main="fragments per RAD site, log10")

# fragments per site for each pool in a scatter plot
	# many sites in pool1 not sequenced, or barely sequenced in pool2
	# these sites may contribute heavily to pool bias. 
plot(log(rowMeans(cot[,keepi & p2]+1),10),log(rowMeans(cot[,keepi & p1]+1),10),ylab="mean coverage per RAD site - pool 1",xlab="mean coverage per RAD site - pool 2")
abline(
	h=c(2,4),
	v=c(0.5,4))

# number non-missing sites for each pool in a scatter plot
	#
	# 
plot(jitter(rowSums(cot[,keepi & p2]>0)),jitter(rowSums(cot[,keepi & p1]>0)),ylab="mean coverage per RAD site - pool 1",xlab="mean coverage per RAD site - pool 2",pch=20,col=rgb(0,0,0,.5))
abline(
	h=c(8.5),
	v=c(5.5))


keepl <- rowSums(cot[,keepi & p2] > 0) > 4 & 
		 rowSums(cot[,keepi & p1] > 0) > 8


# keepl <- log(rowMeans(cot[,keepi & p2]+1),10) > 0.25 & 
# 		 log(rowMeans(cot[,keepi & p2]+1),10) < 4 &
# 		 log(rowMeans(cot[,keepi & p1]+1),10) > 2 &
# 		 log(rowMeans(cot[,keepi & p1]+1),10) < 4

sum(keepl)
rowSums(cot[keepl,keepi]==0) %>% table() %>% plot()
rowSums(cot[keepl,keepi]==0) %>% median()


# replot
plot(log(rowMeans(cot[keepl,keepi & p2]+1),10),log(rowMeans(cot[keepl,keepi & p1]+1),10),ylab="mean coverage per RAD site - pool 1",xlab="mean coverage per RAD site - pool 2")

# total fragments per site for retained individuals and loci
rowSums(cot[keepl,keepi]) %>% log(.,10) %>% hist(.,breaks=50)
rowSums(cot[keepl,keepi]) %>% hist(.,breaks=50)

# get rid of upper tail of coverage. 
	# after this step, 1185 rad sites retained
keepl <- keepl & rowSums(cot[,keepi]) < 100000


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

# mds plot on scaled coverage
md <- (cots[keepl,keepi] + 1) %>% log() %>% t() %>% dist()
mds <- cmdscale(md)
plot(mds,col=factor(substring(rownames(mds),1,5)))

# write list of samples to retain
cat(names(which(keepi)),sep="\n",file="../../../metadata/retained_samples.txt")

# write bed file with intervals to variant call on
bedout <- co[keepl,2:4]
# bed start is zero-indexed
bedout[,2] <- bedout[,2]-1

write.table(bedout,file="../../../metadata/targets.bed",col.names=FALSE,quote=FALSE,row.names=FALSE,sep="\t")


par(mfrow=c(1,3))
ord <- order(unlist(sn[1,]))
plot(unlist(sn[1,])[ord],ylab="raw sequences",col=factor(gsub("_.*","",colnames(sn))[ord]))
abline(h=1e6)
plot(unlist(sn[12,])[ord],col=factor(gsub("_.*","",colnames(sn))[ord]),ylab="reads duplicated")
plot(unlist(sn[26,])[ord],col=factor(gsub("_.*","",colnames(sn))[ord]),ylab="insert size")

plot(unlist(sn[1,])[ord]-unlist(sn[12,])[ord])
