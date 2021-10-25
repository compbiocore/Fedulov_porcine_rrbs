#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -J trim_galore
#SBATCH --mem=32GB
#SBATCH --array=1-18
#SBATCH -e /gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs/rrbs_trimgalore-%A-%a-%J.err
#SBATCH -o /gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs/rrbs_trimgalore-%A-%a-%J.out

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda
conda activate fedulov_rrbs

fastq_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Sequencing_Files"
align_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Alignments"
log_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/Logs"
qc_dir="/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/working/QC"

input=($(ls ${fastq_dir}/*fastq.gz))
bn=$(basename ${input[$((SLURM_ARRAY_TASK_ID -1))]} .fastq.gz)

trim_galore --quality 20 --clip_R1 6 --adapter AGATCGGAAGAGC --output_dir ${fastq_dir} --basename "${bn}_trimgalore" --rrbs ${fastq_dir}/${bn}.fastq.gz

echo ${bn}

