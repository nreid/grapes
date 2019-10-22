#!/bin/bash
#SBATCH --job-name=headcrop
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

module load java
module load Trimmomatic/0.36

# this script trims a single base off the 5' end of R2. 
# there is a spurious T there. 
# it is not removed by stacks. 
# I don't know why. 

# input directory
INDIR=../results/demultiplexed_fastqs

# if alignment dir doesn't exist, make it
mkdir -p ../results/demultiplexed_fastqs/headcrop
OUTDIR=../results/demultiplexed_fastqs/headcrop

# create an array variable containing the file names
FILES=($(ls -1 $INDIR/*.2.fq.gz | grep -v "rem.2"))

INFILE=${FILES[$SLURM_ARRAY_TASK_ID]}
OUTFILE=$(echo $INFILE | sed 's/.*\///')

java -jar $Trimmomatic SE \
-threads 4 \
$INFILE \
$OUTDIR/$OUTFILE \
HEADCROP:1

# mv $OUTDIR/$OUTFILE $INFILE