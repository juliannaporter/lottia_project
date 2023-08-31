#!/bin/bash
#SBATCH --job-name=fasta2bed
#SBATCH -D /ocean/projects/deb200006p/jporter
#SBATCH -o /ocean/projects/deb200006p/jporter/slurm_log/out_%j.txt
#SBATCH -e /ocean/projects/deb200006p/jporter/slurm_log/err_%j.txt
#SBATCH -p RM-shared
#SBATCH -t 48:00:00
#SBATCH --array=0-9

# step one: create bed file from fasta using faidx
fasta=data/genome/Lottia_gigantea.Lotgi1.dna.toplevel.fa
bed=data/genome/Lottia_gigantea.allscaffolds.bed

#pip install pyfaidx
#faidx --transform bed $fasta > $bed

# step two: sort output bedfile
sorted=data/genome/Lottia_gigantea.allscaffolds.sorted.bed

#sort $bed > $sorted

# step three: split sorted bed file using split command
#             (split the same way as delly .excl files)
prefix=data/genome/Lottia_gigantea.splitbed.region

#split -l 450 -d \
# $sorted \
# $prefix

# step four: sort again (just to be consisten w/ delly method)
reg=data/genome/Lottia_gigantea.splitbed.region0${SLURM_ARRAY_TASK_ID}
final=data/genome/Lottia_gigantea.region0${SLURM_ARRAY_TASK_ID}.sorted.bed

#sort $reg > $final

# step five: per manta's instructions, the bed file must be:
#            1) bgzip compressed
#            2) tabix-indexed
module load htslib
bgzip $final
tabix $final.gz
