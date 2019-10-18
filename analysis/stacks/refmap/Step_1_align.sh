#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH -o bwa_%x_%A_%a.out
#SBATCH -e bwa_%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-54]%20


echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

module load samtools/1.7
module load samblaster/0.1.24
module load bwa/0.7.17

INDIR=../results/demultiplexed_fastqs

# create an array variable containing the file names
FILES=($(ls -1 *.R1.fastq))

# get specific file name, assign it to FQ1
	# note that FILES variable is 0-indexed so
	# for convenience we also began the task IDs with 0
FQ1=${FILES[$SLURM_ARRAY_TASK_ID]}
# edit the file name to refer to the mate pair file and assign that name to FQ2
FQ2=$(echo $FQ1 | sed 's/R1/R2/')
# create an output file name
OUT=$(echo $FQ1 | sed 's/.R1.fastq/.sam/')
# write the input filenames to the standard output to check that everything ran according to expectations. 
echo $FQ1 $FQ2

# we won't actually try to align these fake files here but it might look like:

# module load bwa
# bwa mem refgenome.fa $FQ1 $FQ2 > $OUT

echo Files $FQ1 and $FQ2 were aligned by task number $SLURM_ARRAY_TASK_ID on $(date)