#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -n 8
#SBATCH -J trim_galore
#SBATCH --mem=128GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jordan_lawson@brown.edu

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda
conda activate fedulov_rrbs
for sample in `ls /gpfs/data/cbc/fedulov_alexey/porcine_rrbs/Sequencing_Files/*fastq.gz`
do
dir="/gpfs/data/cbc/fedulov_alexey/porcine_rrbs/trimmed"
trim_galore --quality 20 --adapter AGATCGGAAGAGC --output_dir ${dir} --rrbs $sample 
done
