#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-53]%20


echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

module load samtools/1.7
module load samblaster/0.1.24
module load bwa/0.7.17

# input/output directories and supp files
# make sure alignment directory exists

INDIR=../results/demultiplexed_fastqs

# if alignment dir doesn't exist, make it
mkdir -p ../results/aligned_ref
OUTDIR=../results/aligned_ref

# location of genome index
GEN=../../../genome/Vvinifera_145_Genoscope.12X

# create an array variable containing the file names
FILES=($(ls -1 $INDIR/*.1.fq.gz | grep -v "rem.1"))

# get specific file name, assign it to FQ1
	# note that FILES variable is 0-indexed so
	# for convenience we also began the task IDs with 0
FQ1=${FILES[$SLURM_ARRAY_TASK_ID]}
# edit the file name to refer to the mate pair file and assign that name to FQ2
FQ2=$(echo $FQ1 | sed 's/1.fq/2.fq/')
# create an output file name root
SAM=$(echo $FQ1 | sed 's/.1.fq.gz//' | sed 's/.*\///')
# set read group
RG=$(echo \@RG\\tID:$SAM\\tSM:$SAM)

echo $FQ1 $FQ2

# execute the pipe for the son:
bwa mem -t 4 -R $RG $GEN $FQ1 $FQ2 | \
samblaster | \
samtools view -S -h -u - | \
samtools sort -T /scratch/$SAM - >$OUTDIR/${SAM}.bam

echo Files $FQ1 and $FQ2 were aligned by task number $SLURM_ARRAY_TASK_ID on $(date)


