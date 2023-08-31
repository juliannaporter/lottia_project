#!/bin/bash
#SBATCH --job-name=manta
#SBATCH -D /ocean/projects/deb200006p/jporter
#SBATCH -o /ocean/projects/deb200006p/jporter/slurm_log/out_%A_%a.txt
#SBATCH -e /ocean/projects/deb200006p/jporter/slurm_log/err_%A_%a.txt
#SBATCH -p RM
#SBATCH -t 48:00:00
#SBATCH --thread-spec=8
#SBATCH --array=0-9

dir=data/manta/region0${SLURM_ARRAY_TASK_ID}/
ref=data/genome/Lottia_gigantea.Lotgi1.dna.toplevel.fa
reg=data/genome/Lottia_gigantea.region0${SLURM_ARRAY_TASK_ID}.sorted.bed.gz


install/Manta/bin/configManta.py \
 --runDir=$dir \
 --referenceFasta=$ref \
 --callRegions=$reg \
 `cat manta.samples.txt`

data/manta/region0${SLURM_ARRAY_TASK_ID}/runWorkflow.py -j 8
