#!/bin/bash
#SBATCH --job-name=process_radtags
#SBATCH -o process_radtags-%j.out
#SBATCH -e process_radtags-%j.err
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general

# adapted from Madison Caballero
	# specified stacks 2.41. 
	# output set to gzfastq

module load stacks/2.41

#input/output directories, supplementary files
TRIMDIR=../results/trimmed_data/

# make demultiplexed directory if it doesn't exist
mkdir -p ../results/demultiplexed_fastqs
OUTDIR=../results/demultiplexed_fastqs

BARCODES1=../../../metadata/stacks_barcodes_Pool1.txt
BARCODES2=../../../metadata/stacks_barcodes_Pool2.txt


# pool 1
process_radtags \
-1 $TRIMDIR/Pool1_trimmed_1.fastq.gz \
-2 $TRIMDIR/Pool1_trimmed_2.fastq.gz \
-b $BARCODES1 \
-o $OUTDIR/ \
-y gzfastq \
-i gzfastq \
--inline_null \
-e sbfI \
-P \
-c \
-q \
-r \
--barcode_dist_1 1 \
--adapter_1 GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
--adapter_2 AGATCGGAAGAGCGAGAACAA 

# pool 2
process_radtags \
-1 $TRIMDIR/Pool2_trimmed_1.fastq.gz \
-2 $TRIMDIR/Pool2_trimmed_2.fastq.gz \
-b $BARCODES2 \
-o $OUTDIR/ \
-y gzfastq \
-i gzfastq \
--inline_null \
-e sbfI \
-P \
-c \
-q \
-r \
--barcode_dist_1 1 \
--adapter_1 GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
--adapter_2 AGATCGGAAGAGCGAGAACAA 
