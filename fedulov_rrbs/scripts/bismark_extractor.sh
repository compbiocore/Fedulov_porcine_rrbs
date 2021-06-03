#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -N 1 #this says they should all be on the same 
#SBATCH -n 8 #little n in number of tasks
#SBATCH -J methyl_extractor
#SBATCH --mem=20GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jordan_lawson@brown.edu
#SBATCH --array=1-18

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda
conda activate fedulov_rrbs
input=($(ls /gpfs/data/shared/databases/refchef_refs/S_scrofa/primary/alignments/ignore/*.bam)) # using the round brackets indicates that this is a bash array

bismark_methylation_extractor --bedGraph --comprehensive --ignore 3 -s --merge_non_CpG --output /gpfs/data/shared/databases/refchef_refs/S_scrofa/primary/alignments/ignore/ --gzip --multicore 8 --genome_folder /gpfs/data/shared/databases/refchef_refs/S_scrofa/primary/bismark_index/ ${input[$((SLURM_ARRAY_TASK_ID -1))]}  

