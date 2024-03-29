#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -n 32
#SBATCH -J rrbs_genome
#SBATCH --mem=198GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jordan_lawson@brown.edu

source /gpfs/runtime/cbc_conda/bin/activate_cbc_conda
conda activate fedulov_rrbs

mkdir /gpfs/data/shared/databases/refchef_refs/S_scrofa/primary/bismark_index
ln -s /gpfs/data/shared/databases/refchef_refs/S_scrofa/primary/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa /gpfs/data/shared/databases/refchef_refs/S_scrofa/primary/bismark_index/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
bismark_genome_preparation /gpfs/data/shared/databases/refchef_refs/S_scrofa/primary/bismark_index
md5sum /gpfs/data/shared/databases/refchef_refs/S_scrofa/primary/bismark_index/*.* > /gpfs/data/shared/databases/refchef_refs/S_scrofa/primary/bismark_index/final_checksums.md5  
