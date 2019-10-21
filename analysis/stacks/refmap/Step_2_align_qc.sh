#!/bin/bash
#SBATCH --job-name=alignQC
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --qos=general


module load subread/1.6.0
module load samtools/samtools/1.7
module load R/3.6.0
module load bedtools/2.27.1

# This script does some QC on the alignments. 

# make a list of bam files
find ../results/aligned_ref/ -name "*bam" >../results/aligned_ref/bams.list

# make the bam indexes
find ../results/aligned_ref/ -name "*bam" | xargs -I {} samtools index {}

# samtools bam statistics
find ../results/aligned_ref/ -name "*bam" | xargs -I {} samtools stats {}


# find all the restriction sites in the Vvinifera genome assembly
cat ../../../genome/Vvinifera_145_Genoscope.12X.fa | \
Rscript restriction.R > ../../../metadata/vvinifera_sbfi.bed

# buffer restriction sites by 1kb on either side
bedtools slop \
-i ../../../metadata/vvinifera_sbfi.bed \
-g <(cat ../../../genome/Vvinifera_145_Genoscope.12X.fa.fai | cut -f 1-2) \
-b 1000 >../../../metadata/vvinifera_int.bed

# run featureCounts to get a sense of the number of fragments on RAD sites

# first generate a SAF file for featureCounts:
	# use dummy gene IDs
echo "GeneID	Chr	Start	End	Strand" >../../../metadata/vvinifera_int.saf
cat ../../../metadata/vvinifera_int.bed | \
awk '{OFS="\t"}{s=$2+1}{print NR,$1,s,$3,"+"}' >>../../../metadata/vvinifera_int.saf

# run featurecounts
mkdir -p ../results/align_stats
featureCounts \
-a ../../../metadata/vvinifera_int.saf \
-o ../results/align_stats/vvinifer_counts.txt \
-Q 30 \
-F SAF \
--primary \
-p \
-T 12 \
$(cat ../results/aligned_ref/bams.list | tr "\n" " ")




