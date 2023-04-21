#!/bin/bash -l
#SBATCH -A naiss2023-22-266
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 39:00:00
#SBATCH -J count
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL

cd /crex/proj/snic2021-23-14/Xue/scRNA-seq/fastq/GSE123013
module load bioinfo-tools cellranger/7.0.1
cellranger count --id=GSE123013 \
--transcriptome=/crex/proj/snic2021-23-14/Xue/Ath \
--fastqs=/crex/proj/snic2021-23-14/Xue/scRNA-seq/fastq/GSE123013 \
--sample=WT-WERGFP,WT-WERGFP_2,WT-WERGFP_3 \
--localcores=8 \
--localmem=64