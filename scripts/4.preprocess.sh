#!/bin/bash
#SBATCH --job-name=bcftools
#SBATCH -D /ocean/projects/deb200006p/jporter
#SBATCH -o /ocean/projects/deb200006p/jporter/slurm_log/out_%A_%a.txt
#SBATCH -e /ocean/projects/deb200006p/jporter/slurm_log/err_%A_%a.txt
#SBATCH -p RM-shared
#SBATCH -t 48:00:00
#SBATCH --array=0-9

ref=data/genome/Lottia_gigantea.Lotgi1.dna.toplevel.fa

init_d=data/delly/ALL.delly.call.region0${SLURM_ARRAY_TASK_ID}.vcf
init_m=data/manta/region0${SLURM_ARRAY_TASK_ID}/results/variants/diploidSV.vcf.gz

# step 1: process files for merging with bcftools merge
# delly
norm_d=data/delly/ALL.delly.region0${SLURM_ARRAY_TASK_ID}.norm.vcf
tbi_d=data/delly/ALL.delly.region0${SLURM_ARRAY_TASK_ID}.norm.vcf.gz

module load bcftools
bcftools norm --check-ref ws \
 -Ov -o $norm_d \
 --fasta-ref $ref \
 $init_d

module load htslib
bgzip $norm_d
tabix $tbi_d

# manta
tbi_m=data/manta/ALL.manta.region0${SLURM_ARRAY_TASK_ID}.vcf.gz

#mv $init_m $tbi_m
#tabix $tbi_m

