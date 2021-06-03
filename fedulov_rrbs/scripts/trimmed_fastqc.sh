#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -n 8
#SBATCH -J trimmed_fastqc
#SBATCH --mem=16GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jordan_lawson@brown.edu

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda
conda activate fedulov_rrbs
for sample in `ls /gpfs/data/shared/databases/refchef_refs/S_scrofa/primary/trimmed/*_trimmed.fq.gz`
do
trim_qc_dir="/gpfs/data/shared/databases/refchef_refs/S_scrofa/primary"
fastqc -o ${trim_qc_dir}/trimmed_qc $sample
done
