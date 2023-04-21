#!/bin/bash -l
#SBATCH -A naiss2023-22-266
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 27:00:00
#SBATCH -J index
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL
module load bioinfo-tools cellranger/7.0.1
cd /crex/proj/snic2021-23-14/Xue
cellranger mkref \
--genome=Ath \
--fasta=Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
--genes=Ath.gtf