#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -J bismark_align_update
#SBATCH --mem=160GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jordan_lawson@brown.edu
#SBATCH --array=1-18
#SBATCH -e /gpfs/data/cbc/fedulov_alexey/porcine_rrbs/logs/bismark_align_%a_%A_%j.err
#SBATCH -o /gpfs/data/cbc/fedulov_alexey/porcine_rrbs/logs/bismark_align_%a_%A_%j.out

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda
conda activate fedulov_rrbs
input=($(ls /gpfs/data/cbc/fedulov_alexey/porcine_rrbs/trimmomatic_update/*_tr.fq.gz)) # using the round brackets indicates that this is a bash array
bismark -o /gpfs/data/cbc/fedulov_alexey/porcine_rrbs/alignments/update --bowtie2 --genome /gpfs/data/shared/databases/refchef_refs/S_scrofa/primary/bismark_index/ --un --pbat ${input[$((SLURM_ARRAY_TASK_ID -1))]}
