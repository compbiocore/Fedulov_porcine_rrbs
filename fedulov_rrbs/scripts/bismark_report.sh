#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -N 1 #this says they should all be on the same 
#SBATCH -n 8 #little n in number of tasks
#SBATCH -J bismark_summary
#SBATCH --mem=10GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jordan_lawson@brown.edu
#SBATCH --array=1-18

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda
conda activate fedulov_rrbs
input=($(ls /gpfs/data/cbc/fedulov_alexey/porcine_rrbs/alignments/ignore/*_bismark_bt2_SE_report.txt)) # using the round brackets indicates that this is a bash array
bismark2report ${input[$((SLURM_ARRAY_TASK_ID -1))]}  
