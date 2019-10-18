#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH -o trimmomatic-%j.out
#SBATCH -e trimmomatic-%j.err
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general

# adapted from Madison Caballero

module load java
module load Trimmomatic/0.36

#input/output directories, supplementary files
POOL1=/labs/Wegrzyn/Grapes/data/raw_data/1
POOL2=/labs/Wegrzyn/Grapes/data/raw_data/2
ADAPT=../../../metadata/adapters.fa

mkdir -p ../results/trimmed_data
OUTDIR=../results/trimmed_data


java -jar $Trimmomatic PE \
$POOL1/Pool1_1.fastq.gz \
$POOL1/Pool1_2.fastq.gz \
Pool1_trimmed_1.fastq.gz \
unpaired_Pool1_trimmed_1.fastq.gz \
Pool1_trimmed_2.fastq.gz \
unpaired_Pool1_trimmed_2.fastq.gz \
ILLUMINACLIP:$ADAPT:3:30:10 \
SLIDINGWINDOW:3:10 \
MINLEN:50 


java -jar $Trimmomatic PE \
$POOL2/Pool2_1.fastq.gz \
$POOL2/Pool2_2.fastq.gz \
Pool2_trimmed_1.fastq.gz \
unpaired_Pool2_trimmed_1.fastq.gz \
Pool2_trimmed_2.fastq.gz \
unpaired_Pool2_trimmed_2.fastq.gz \
ILLUMINACLIP:$ADAPT:3:30:10 \
SLIDINGWINDOW:3:10 \
MINLEN:50

