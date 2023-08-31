#!/bin/bash
#SBATCH --job-name=bcftools
#SBATCH -D /ocean/projects/deb200006p/jporter
#SBATCH -o /ocean/projects/deb200006p/jporter/slurm_log/out_%j.txt
#SBATCH -e /ocean/projects/deb200006p/jporter/slurm_log/err_%j.txt
#SBATCH -p RM
#SBATCH -t 48:00:00

module load bcftools

#bcftools +split \
# -o data/delly/delly.split.wholegenome.sample \
# data/delly/All.delly.whole.genome.rename.vcf

bcftools +split \
 -o data/manta/manta.split.wholegenome.sample \
 data/manta/All.manta.whole.genome.rename.vcf
