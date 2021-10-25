#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -J fastqc_trimmed
#SBATCH --mem=32GB
#SBATCH --array=1-18
#SBATCH -e /gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs/rrbs_fastqc_trimmed-%A-%a-%J.err
#SBATCH -o /gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs/rrbs_fastqc_trimmed-%A-%a-%J.out

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda

fastq_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Sequencing_Files"
align_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Alignments"
log_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs"
qc_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/QC"

input=($(ls ${fastq_dir}/*_trimmomatic.fastq.gz))
bn=$(basename ${input[$((SLURM_ARRAY_TASK_ID -1))]} _trimmomatic.fastq.gz)

fastqc -o ${qc_dir}/fastqc ${fastq_dir}/${bn}_trimmomatic.fastq.gz
