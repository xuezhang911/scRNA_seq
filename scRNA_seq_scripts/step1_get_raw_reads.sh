#!/bin/bash -l
#SBATCH -A naiss2023-22-266
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 29:00:00
#SBATCH -J getreads
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL

cd /crex/proj/snic2021-23-14/Xue/scRNA-seq/fastq/GSE123013
module load bioinfo-tools sratools/3.0.0
for i in $(seq 100 103);
do fastq-dump  --split-files --gzip  SRR8257$i
done

