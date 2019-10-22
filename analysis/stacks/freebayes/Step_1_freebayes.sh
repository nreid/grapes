#!/bin/bash 
#SBATCH --job-name=freebayes
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 7
#SBATCH --mem=30G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

# load required software

module load bedtools
module load bamtools
module load htslib
module load freebayes


# make a directory for results if it doesn't exist
mkdir -p ../results/freebayes 
OUTDIR=../results/freebayes

# make a list of bam files
	# even though we've filtered individuals, try to call variants on all of them. 
find ../results/aligned_ref/ -name "*bam" >../results/aligned_ref/bams.list

# set a variable for the reference genome location
GEN=../../../genome/Vvinifera_145_Genoscope.12X.fa

# set a variable for the rad sites targeted
TARGETS=../../../metadata/targets.bed

# set a variable for the bam list
BAMLIST=../results/aligned_ref/bams.list

# note that bamtools region specification uses ".." instead of "-"
bamtools merge -list $BAMLIST | \
bamtools filter -in stdin -mapQuality ">30" -isProperPair true | \
bedtools intersect -a stdin -b $TARGETS -nonamecheck | \
freebayes -f $GEN --stdin | \
bgzip -c >$OUTDIR/vvinifera_fb.vcf.gz

tabix -p vcf $OUTDIR/vvinifera_fb.vcf.gz

date
