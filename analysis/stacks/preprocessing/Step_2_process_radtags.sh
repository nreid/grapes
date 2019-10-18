#!/bin/bash
#SBATCH --job-name=process_radtags
#SBATCH -o process_radtags-%j.o
#SBATCH -e process_radtags-%j.e
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general

# adapted from Madison Caballero

module load stacks/2.41

# input/output directories
TRIMDIR=../trimmed_data/
OUTDIR=demultiplex_output

# pool 1
process_radtags \
-1 $TRIMDIR/Pool1_trimmed_1.fastq.gz \
-2 $TRIMDIR/Pool1_trimmed_2.fastq.gz \
-b barcodes_Pool1 \
-o $OUTDIR/ \
-y fastq \
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
-b barcodes_Pool2 \
-o $OUTDIR/ \
-y fastq \
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
