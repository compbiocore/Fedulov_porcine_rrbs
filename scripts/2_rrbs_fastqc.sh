#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -J rrbs_fastqc
#SBATCH --mem=16GB
#SBATCH --mail-type=ALL
#SBATCH --array=1-18
#SBATCH -e /gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs/rrbs_fastqc-%A-%a-%J.err
#SBATCH -o /gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs/rrbs_fastqc-%A-%a-%J.out

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda
conda activate fedulov_rrbs

fastq_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Sequencing_Files"
align_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Alignments"
log_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs"
qc_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/QC"

input=($(ls ${fastq_dir}/*_R1_001.fastq.gz))
bn=$(basename ${input[$((SLURM_ARRAY_TASK_ID -1))]} _R1_001.fastq.gz)

fastqc -o ${qc_dir}/fastqc ${fastq_dir}/${bn}_R1_001.fastq.gz
