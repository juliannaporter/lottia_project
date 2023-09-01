#!/bin/bash
#SBATCH --job-name=split_regions
#SBATCH -D /ocean/projects/deb200006p/jporter
#SBATCH -o /ocean/projects/deb200006p/jporter/slurm_log/out_%j.txt
#SBATCH -e /ocean/projects/deb200006p/jporter/slurm_log/err_%j.txt
#SBATCH -p RM-shared
#SBATCH -t 48:00:00
#SBATCH --array=0-9

#split -l 450 -d \
# data/genome/Lottia_gigantea.scafname.sorted \
# data/genome/Lottia_gigantea.splitscaf.region

sort data/genome/Lottia_gigantea.splitscaf.region0${SLURM_ARRAY_TASK_ID} > \
 data/genome/Lottia_gigantea.splitscaf.region0${SLURM_ARRAY_TASK_ID}.sorted

comm -23 data/genome/Lottia_gigantea.scafname.sorted \
 data/genome/Lottia_gigantea.splitscaf.region0${SLURM_ARRAY_TASK_ID}.sorted \
 > data/genome/Lottia_gigantea.region0${SLURM_ARRAY_TASK_ID}.excl

comm -23 data/genome/
