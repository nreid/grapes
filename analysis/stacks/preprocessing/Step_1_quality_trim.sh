#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH -o trimmomatic-%j.out
#SBATCH -e trimmomatic-%j.err
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --partition=general
#SBATCH --qos=general

# adapted from Madison Caballero
	# reduced quality trimming stringency to meanQ < 10

module load java
module load Trimmomatic/0.36

#input/output directories, supplementary files
POOL1=/labs/Wegrzyn/Grapes/data/raw_data/1
POOL2=/labs/Wegrzyn/Grapes/data/raw_data/2
ADAPT=../../../metadata/adapters.fa

# make trimmed directory if it doesn't exist
mkdir -p ../results/trimmed_data
OUTDIR=../results/trimmed_data

# pool 1
java -jar $Trimmomatic PE \
-threads 20 \
$POOL1/Pool1_R1_.fastq.gz \
$POOL1/Pool1_R2_.fastq.gz \
$OUTDIR/Pool1_trimmed_1.fastq.gz \
$OUTDIR/unpaired_Pool1_trimmed_1.fastq.gz \
$OUTDIR/Pool1_trimmed_2.fastq.gz \
$OUTDIR/unpaired_Pool1_trimmed_2.fastq.gz \
ILLUMINACLIP:$ADAPT:3:30:10 \
SLIDINGWINDOW:3:10 \
MINLEN:50 

# pool 2
java -jar $Trimmomatic PE \
-threads 20 \
$POOL2/Pool2_R1_.fastq.gz \
$POOL2/Pool2_R2_.fastq.gz \
$OUTDIR/Pool2_trimmed_1.fastq.gz \
$OUTDIR/unpaired_Pool2_trimmed_1.fastq.gz \
$OUTDIR/Pool2_trimmed_2.fastq.gz \
$OUTDIR/unpaired_Pool2_trimmed_2.fastq.gz \
ILLUMINACLIP:$ADAPT:3:30:10 \
SLIDINGWINDOW:3:10 \
MINLEN:50

