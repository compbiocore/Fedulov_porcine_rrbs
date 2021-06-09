#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -n 32
#SBATCH -J trimmomatic
#SBATCH --mem=198GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jordan_lawson@brown.edu

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda

for sample in `ls /gpfs/data/cbc/fedulov_alexey/porcine_rrbs/trimmed/*_trimmed.fq.gz`
do
    dir="/gpfs/data/cbc/fedulov_alexey/porcine_rrbs/trimmomatic"
    base=$(basename $sample "_trimmed.fq.gz")
    trimmomatic SE  -threads 8 -trimlog ${dir}/${base}_SE.log $sample ${dir}/${base}_tr.fq.gz ILLUMINACLIP:/gpfs/data/cbc/cbc_conda_v1/envs/cbc_conda/opt/trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:5:6:true SLIDINGWINDOW:4:20 MINLEN:35
done 


