#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -J trimmomatic
#SBATCH --mem=32GB
#SBATCH --array=1-18
#SBATCH -e /gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs/rrbs_trimmomatic-%A-%a-%J.err
#SBATCH -o /gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs/rrbs_trimmomatic-%A-%a-%J.out

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda

fastq_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Sequencing_Files"
align_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Alignments"
log_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs"
qc_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/QC"

input=($(ls ${fastq_dir}/*_trimgalore_trimmed.fq.gz))
bn=$(basename ${input[$((SLURM_ARRAY_TASK_ID -1))]} _trimgalore_trimmed.fq.gz)

trimmomatic SE -threads 8 -trimlog ${log_dir}/${bn}_trimmomatic.log ${fastq_dir}/${bn}_trimgalore_trimmed.fq.gz ${fastq_dir}/${bn}_trimmomatic.fastq.gz ILLUMINACLIP:/gpfs/data/cbc/cbc_conda_v1/envs/cbc_conda/opt/trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:5:6:true SLIDINGWINDOW:4:20 MINLEN:35

