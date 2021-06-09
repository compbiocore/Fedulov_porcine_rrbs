# Fedulov RRBS porcine data analysis pipeline 

## Data

All of the sequencing files for this project can be found in `/gpfs/data/cbc/fedulov_alexey/porcine_rrbs/Sequencing_Files`. 

# Converting Reference Genome

First, I made the Bismark genome index using: 

```{bash}
sbatch rrbs_genome.sh 
```

## Initial QC

Secondly, I ran a few QC metrics using the `rrbs_fastqc.sh` script, which will run fastqc on the fastq sequencing files: 

```{bash}
sbatch rrbs_fastqc.sh
```
## Trimming via Trim Galore

Next, I used trim galore in the RRBS mode to trim the reads:

```{bash}
sbatch trim_galore.sh
```   
Then I ran the reads through fastqc:

```{bash}
sbatch trimmed_fastqc.sh
```  

## Trimmomatic 

I then ran trimmomatic on the trimmed reads to remove TruSeq adapters: 

```{bash}
sbatch trimmomatic.sh
```  

And once again ran the reads through fastqc: 

```{bash}
sbatch retrim_fastqc.sh
``` 

## Bismark aligning

Map the reads using the `--non_directional` flag: 

```{bash}
sbatch bismark_align.sh
```

Subsequently, run Bismark extractor (note: the code in this step still needs to be edited based on results from previous step):

```{bash}
bismark_extractor.sh
```

Make some reports:

```{bash}
sbatch bismark_report.sh
```
```{bash}
sbatch bam2nuc.sh
```

Note: Add steps and edit as necessary





