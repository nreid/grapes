# grapes

repository contains two projects using armenian grape data
- a rad-seq tutorial 
	- stacks
	- ipyrad
- a test case for nextflow

___

raw data located in /labs/Wegrzyn/Grapes/data/raw_data

54 samples in 2 pools of 27 individuals. 

___

scripts adapted from Madison Caballero

___

analysis indicates this data set has some problems:

1. median pcr duplication rate 92%
2. spurious T beginning read 2. 
3. large batch effects
	- insert sizes differ between pools
	- bias in coverage of rad sites between pools
4. highly unequal pooling of individuals
5. high levels of missing data

salvage strategy:

1. filter individuals and rad sites with high stringency on missing data rate
2. haploidize genotypes
3. set high minor allele frequency threshold