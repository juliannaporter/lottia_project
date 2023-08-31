#!/bin/bash
#SBATCH --job-name=rename
#SBATCH -D /ocean/projects/deb200006p/jporter
#SBATCH -o /ocean/projects/deb200006p/jporter/slurm_log/out_%j.txt
#SBATCH -e /ocean/projects/deb200006p/jporter/slurm_log/err_%j.txt
#SBATCH -p RM-shared
#SBATCH -t 48:00:00

module load bcftools

bcftools reheader \
 --samples data/delly/sample.rename \
 --output data/delly/All.delly.whole.genome.rename.vcf \
 data/delly/All.delly.whole.genome.vcf

bcftools reheader \
 --samples data/delly/sample.rename \
 --output data/manta/All.manta.whole.genome.rename.vcf \
 data/manta/All.manta.whole.genome.vcf
