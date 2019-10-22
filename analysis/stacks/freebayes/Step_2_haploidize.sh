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


module load htslib/1.7
module load R/3.6.0

zcat vvinifera_fb.vcf.gz | \
Rscript haploidize_stream.R | \
bgzip >vvinifera_fb_hap.vcf.gz

tabix -p vcf vvinifera_fb_hap.vcf.gz
