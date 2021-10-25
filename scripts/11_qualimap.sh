#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -N 1 
#SBATCH -n 8 
#SBATCH -J qualimap
#SBATCH --mem=8GB
#SBATCH --array=1-18
#SBATCH -e /gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs/qualimap-%A-%a-%J.err
#SBATCH -o /gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs/qualimap-%A-%a-%J.out

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda
conda activate fedulov_rrbs

fastq_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Sequencing_Files"
align_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Alignments"
log_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs"
qc_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/QC"

input=($(ls ${align_dir}/*_trimmomatic_bismark_bt2.sorted.bam)) 
bn=$(basename ${input[$((SLURM_ARRAY_TASK_ID -1))]} _trimmomatic_bismark_bt2.sorted.bam)

qualimap bamqc -bam ${align_dir}/${bn}_trimmomatic_bismark_bt2.sorted.bam -outdir ${qc_dir}/qualimap

