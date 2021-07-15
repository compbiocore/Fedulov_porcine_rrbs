#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -n 32
#SBATCH -J qualimap
#SBATCH --mem=198GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jordan_lawson@brown.edu

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda

for sample in `ls /gpfs/data/cbc/fedulov_alexey/porcine_rrbs/alignments/update/*_tr_bismark_bt2.bam`
do
    dir="/gpfs/data/cbc/fedulov_alexey/porcine_rrbs/qualimap"
    base=$(basename $sample "_tr_bismark_bt2.bam")
    samtools sort -o ${dir}/${base}.srtd.bam $sample
    qualimap bamqc -bam ${dir}/${base}.srtd.bam -outdir ${dir}/qc_${base}
done 
