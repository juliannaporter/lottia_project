#!/bin/bash
#SBATCH --job-name=delly
#SBATCH -D /ocean/projects/deb200006p/jporter
#SBATCH -o /ocean/projects/deb200006p/jporter/slurm_log/out_%j.txt
#SBATCH -e /ocean/projects/deb200006p/jporter/slurm_log/err_%j.txt
#SBATCH -p RM
#SBATCH -t 48:00:00
##SBATCH --array=0-9
#SBATCH --thread-spec=8

make PARALLEL=1 -B install/delly/delly
export OMP_NUM_THREADS=8

#dir=../enielsen/LGwork/08_realigned
#dir=data/bam
#CNV=$(ls $dir/*.bam.realigned.bam)
#bam="ls $dir/*bam.realigned.bam"

ref=data/genome/Lottia_gigantea.Lotgi1.dna.toplevel.fa

out1=data/delly/ALL.delly.call.whole.genome.vcf
out2=data/delly/ALL.delly.call.region0${SLURM_ARRAY_TASK_ID}.vcf

region=data/genome/Lottia_gigantea.region0${SLURM_ARRAY_TASK_ID}.excl

install/delly/delly call \
 -g $ref \
 -t DEL,DUP,INS \
 -q 20 \
 `echo data/bam/*.bam` > $out1



#install/delly/delly call \
# -g $ref \
# -t DEL,DUP,INS \
# -q 20 \
# -x $region \
# `echo data/bam/*.bam` > $out2

