# Fedulov RRBS porcine data analysis pipeline 

This file contains the steps for re-producing the porcine rrbs analysis, outlining the order in which the scripts were ran and what was performed for downstream analyses. All scripts for the analyses can be found in the scripts folder of this repo. 

## Data

All of the sequencing files for this project can be found in `/gpfs/data/cbc/fedulov_alexey/porcine_rrbs/Sequencing_Files`. 

# Converting Reference Genome

First, we made the Bismark genome index using: 

```{bash}
sbatch rrbs_genome.sh 
```

## Initial QC

Secondly, we ran a few QC metrics using the `rrbs_fastqc.sh` script, which will run fastqc on the fastq sequencing files (i.e., the raw reads): 

```{bash}
sbatch rrbs_fastqc.sh
```
## Trimming via Trim Galore

Next, we used trim galore in the RRBS mode to trim the reads:

```{bash}
sbatch trim_galore_5prime.sh
```   
Then we ran these trimmed reads through fastqc using:

```{bash}
sbatch trimmed_fastqc_update.sh
```  

## Trimmomatic 

Subsequently, we performed a second round of trimming to remove TruSeq adapters using the following: 

```{bash}
sbatch trimmomatic_update.sh
```  

And once again ran the trimmed reads through fastqc: 

```{bash}
sbatch retrim_fastqc_update.sh
``` 

## Bismark aligning and report generation 

We then mapped the reads using the `--non_directional` flag with the following script: 

```{bash}
sbatch bismark_align_update.sh
```

Subsequently, we performed methylation extraction using Bismark extractor:

```{bash}
bismark_extractor.sh
```

We then make some reports with the following scripts:

```{bash}
sbatch bismark_report.sh
```
```{bash}
sbatch bam2nuc.sh
```

# Downstream R analysis 

After this point, we used the bam coverage files to perform analyses in R. The code for the R analysis can be found in the `R_rrbs_updated analysis.R` file. For our downstream analysis, genes were filtered such that only CpGs with at least 5 counts (methylated and unmethylated) in 3 samples were included for downstream analysis. Additionally, CpGs that were never methylated or always methylated were filtered out, as these provide no information about differential methylation. The resulting sample size after filtering was 84,999.

We then used the glmFit function to fit a negative binomial generalized log-linear model. The experimental design matrix was constructed using modelMatrixMeth with a factorial experimental design (~0 + group), where group was a factor variable with levels comprised of each combination of treatment/diet/tissue. We dropped the intercept from our model to re-parameterize it as a means model. Subsequently, the glmLRT function was used to find differentially methylated loci for comparisons of interest, which were made by constructing contrast vectors. Individual CpG sites were considered differentially methylated if the nominal p-value was < 0.01. Results were filtered based on this nominal p-value threshold.





