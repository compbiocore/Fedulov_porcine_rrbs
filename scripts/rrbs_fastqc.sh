#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -n 32
#SBATCH -J rrbs_fastqc
#SBATCH --mem=198GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jordan_lawson@brown.edu

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda
conda activate fedulov_rrbs
for sample in `ls /gpfs/data/cbc/fedulov_alexey/porcine_rrbs/Sequencing_Files/*fastq.gz`
do
align_dir="/gpfs/data/cbc/fedulov_alexey/porcine_rrbs" 
fastqc -o ${align_dir}/fastqc $sample
done
