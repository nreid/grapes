# Stacks

54 samples in 2 pools. 

file 'metadata.txt' was written by hand in excel originally. 
I edited it to make it friendlier, but still must be read in R as:
`read.table("metadata.txt",header=TRUE,sep="\t",comment.char="",quote="",stringsAsFactors=FALSE)`

## Process pooled data

1. Trimmomatic
2. stacks: process_radtags

## Reference mapping approach

Must have bwa indexed genome in repository root dir `genome`. 

1. Align to reference genome
2. Explore results of mapping

### Variant calling

1. Inserted a freebayes variant call in here. Not stacks, but...