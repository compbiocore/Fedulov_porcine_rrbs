# Fedulov RRBS porcine data analysis pipeline 

This file contains the steps for re-producing the porcine rrbs analysis, outlining the order in which the scripts were ran and what was performed for downstream analyses. All scripts for the analyses can be found in the scripts folder.     


### Upstream analyses in Conda

The upstream analysis portion of this analysis (everything up to and including the creation of the Bismark coverage files) was ran on Oscar in a conda environment. The conda environment can be re-created using the `dat/porcine_rrbs_conda.yml` file. 

```{bash}
conda env create -f porcine_rrbs_conda.yaml
```

And then activating the conda environment:
```{bash}
conda activate fedulov_rrbs
```

### Downstream analyses in Docker

The downstream analyses were all run in a Docker container, which can be accessed by pulling the docker image:

```{bash}
docker pull compbiocore/woodsc:latest
```

And then running the container:

```{bash}
docker run --rm -p 8787:8787 -e USER=rstudio -e PASSWORD=yourpassword --volume ${PWD}:/home/rstudio rrbs:latest 
```

Opening a browser and navigating to `localhost:8787` and entering the USER and PASSWORD as indicated above will start up RStudio with all packages installed necessary for the downstream analysis.

## Data Organization

The raw data and intermediate analysis files are on Oscar in the `/gpfs/data/cbc/fedulov_alexey/porcine_rrbs` directory.    


All of the raw sequencing files for this project can be found on Oscar in `/gpfs/data/cbc/fedulov_alexey/porcine_rrbs/Sequencing_Files`. 


## Converting Reference Genome

First, we made the Bismark genome index using the `bismark_genome_preparation` function in Bismark: 

```{bash}
sbatch rrbs_genome.sh 
```

## Initial QC

Secondly, we ran a few QC metrics on the fastq sequencing files (i.e., the raw reads): 

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

## Downstream R analysis 

After this point, we used the bam coverage files to perform analyses in R. The code for the R analysis can be found in the `R_rrbs_updated analysis.R` file. For our downstream analysis, the data was filtered so that each CpG loci with at least 5 counts (methylated and unmethylated) in 3 samples were included for downstream analysis and removed any CpGs that were never or always methylated. We also removed any non-chromosomal scaffolds from our analysis. 

The biomaRt package in R was used to find all transcripton start sites (TSS) in the S. scrofa genome. The `DistanceToNearest` function was used to determine each CpGs distance to its nearest TSS and loci that were more than 2 kb up or down stream from a TSS were removed from analysis. CpG that were within 2KB of a TSS were considered 'promoters' and reads were summed across each promoter. Further downstream analyses used these summed counts per promoter regions as their unit of measurement.

We then used the `glmFit` function to fit a negative binomial generalized log-linear model. The experimental design matrix was constructed using `modelMatrixMeth` with a factorial experimental design (~0 + group), where group was a factor variable with levels comprised of each combination of treatment/diet/tissue. We dropped the intercept from our model to re-parameterize it as a means model. Subsequently, the `glmLRT` function was used to find differentially methylated loci for comparisons of interest, which were made by using the `makeContrasts` function. Individual CpG sites were considered differentially methylated if the nominal p-value was < 0.01. Results were filtered based on this nominal p-value threshold.

The following contrasts were ran:

`SMVisch_vs_SMVnormal` uses SMV normal as the reference group and contrasts it with the SMV ischemic group.





