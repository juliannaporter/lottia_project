#!/bin/bash
#SBATCH --job-name=merge
#SBATCH -D /ocean/projects/deb200006p/jporter
#SBATCH -o /ocean/projects/deb200006p/jporter/slurm_log/out_%j.txt
#SBATCH -e /ocean/projects/deb200006p/jporter/slurm_log/err_%j.txt
#SBATCH -p RM-shared
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=8

module load bcftools

#for caller in manta delly
for caller in delly
do
 output=data/$caller/All.$caller.whole.genome.vcf
 bcftools concat \
 -o $output -Ov \
 `echo data/$caller/ALL.$caller.region*.gz`
done
