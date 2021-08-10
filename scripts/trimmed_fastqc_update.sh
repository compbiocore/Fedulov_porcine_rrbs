#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -n 8
#SBATCH -J trimmed_fastqc_update
#SBATCH --mem=16GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jordan_lawson@brown.edu

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda
conda activate fedulov_rrbs
for sample in `ls /gpfs/data/cbc/fedulov_alexey/porcine_rrbs/trim_galore/*_trimmed.fq.gz`
do
trim_qc_dir="/gpfs/data/cbc/fedulov_alexey/porcine_rrbs"
fastqc -o ${trim_qc_dir}/trim_galore_qc $sample
done
