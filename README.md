# Fedulov Porcine RRBS

This repository tracks the analysis pipeline for the RRBS analysis of porcine data for the Fedulov lab. All analyses were ran on Oscar in the folder `/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs`.

# Organization overview
- The `scripts` folder contains all of the scripts used to run the pipeline (from raw reads to Bismark coverage files). These scripts use a Conda environment on Oscar.
- The `metadata` folder contains the sample manifest (`targets.txt`, the yml file to create the computing environment for upstream analyses (`porcine_rrbs_conda.yml`), the yaml used to document the commands to create the Bismark references (`new.yaml`), the `Dockerfile` to create the environment for downstream analyses and the `packages.R` file to install R packages in the Docker container.
- The `notebooks` folder contains the notebooks for downstream analyses (from Bismark coverage files to differentially methylated loci). These notebooks were run on Oscar inside a Singularity container.
- The `working` directory on Oscar contains folders and files that are not pushed to git. This includes folders from the sequencing centers initial attempt at the analyses (`epigentek`) and the raw data (`Sequencing_files`). It also includes our read QC (`fastqc`), trimmed reads (`trim_galore` and `trimmomatic`), trimmed read QC (`trim_galore_qc` and `trimmomatic_qc`), Bismark alignments and reports (`alignments`), alignment qc (`qualimap`), and log files (`logs`).

# Setting up the environment

Running the contents of `scripts` directory on Oscar shouldn't require any extra set up, as they reference the proper conda environments when they are executed. You can also recreate the computing environment using ` metadata/porcine_rrbs_conda.yml`

To run the contents of `notebooks`, you can use Singularity on Oscar (or Docker on your local machine). To run the Docker container on your local machine, clone this repository and `cd` into it and run:

```{bash}
docker pull compbiocore/rrbs:latest .
```

To run the container and mount the present working directory to `/home/rstudio`:

```{bash
docker run --rm -p 8787:8787 -e USER=rstudio -e PASSWORD=passwordhere --volume ${PWD}:/home/rstudio rrbs:latest
```

Then navigate to localhost:8787 in a browser.


To run as Singularity on Oscar:

```{bash}
cd /gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/metadata
mkdir -p run var-lib-rstudio-server
printf 'provider=sqlite\ndirectory=/var/lib/rstudio-server\n' > database.conf
singularity pull docker://compbiocore/rrbs:latest
```

Then log into oscar over VNC client, then open terminal and run:

```{bash}
cd /gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs/metadata
export SINGULARITY_BINDPATH="/gpfs/data/cbc"
singularity exec --bind run:/run,var-lib-rstudio-server:/var/lib/rstudio-server,database.conf:/etc/rstudio/database.conf rrbs_latest.sif rserver --www-address=127.0.0.1
```

Then open another terminal window and run:
```{bash}
module load firefox/87.0
firefox
```

Then navigate to localhost:8787

Once R is open, you'll need to tell R to use the libraries installed in the singularity container (rather than looking in your home directory).

```{r}
.libPaths(c('/usr/local/lib/R/site-library', '/usr/local/lib/R/library'))
```

# Experimental Design

Study design: Pigs received either normal diet or high-fat diet; furthermore, in either diet group we have normal heart tissue vs. surgically induced ischemic heart tissue, 
and in either diet group the ischemia was treated with microvesicles injection.

The following 6 groups (n=3 each) of porcine myocardium DNA are therefore presented: 
1. “SMV normal” = sham surgery with normal diet (normal heart tissue, negative control);
2. “SMV ischemic” = sham surgery with normal diet (ischemic heart tissue); 
3. “MVM” = myocardial injection of microvesicles in a normal diet model; 
4. “HVM” = myocardial injection of microvesicles in a high-fat diet model; 
5. “HSMV normal” = myocardial injection of vehicle control in a high-fat diet - normal heart tissue; 
6. “HSMV ischemic” = myocardial injection of vehicle control in a high-fat diet – ischemic heart tissue.

Goal: we want to determine the loci differentially methylated in all pairwise comparisons between the groups to answer these experimental questions:

1). What was the effect of ischemia vs. normal heart tissue on a normal diet in absence of treatments? (SMV normal vs. SMV ischemic).

2). What was the effect of ischemia vs. normal heart tissue on a high-fat diet in absence of treatments (HSMV normal vs. HSMV ischemic) and whether it was different from that effect in the normal diet? (this resultant list vs. list from #1) 

3). What was the effect of microvesicle treatment (vs. ischemia) in normal diet and in high-fat diet? (MVM vs. SMV ischemic; and HVM vs. HSMV ischemic)

Notes: The differentially methylated loci detection should be done at reasonable stringency to a) characterize the extent of changes by each experimental factor, and b) to inform downstream pathway analysis.

The differentially methylated loci need to be mapped to the nearest transcript (with distance to TSS) and presented in the form of data tables convenient for heatmaps (i.e. with methylation percent values) and for downstream pathway analysis (i.e. with gene names and IDs).

# Analysis Overview
Initial QC of raw reads was run using FASTQC (0.11.9) <sup>1</sup>. Reads were trimmed using Trim Galore (0.6.5) <sup>2</sup> with the `--rrbs` option to remove any filled-in cytosine positions remaining from the library preparation, `--adapter AGATCGGAAGAGC` to remove adapters, `--clip_R1 6` to remove primers left from the PBAT library preparation, and `--quality 20` to remove low-quality ends from the reads. An additional trimming step was performed in Trimmomatic (0.36) to remove any additional adapters (`ILLUMINACLIP:TruSeq3-SE.fa:2:30:5:6:true`), low-quality bases (`SLIDINGWINDOW:4:20`), or truncated reads (`MINLEN:35`). Bismark (0.22.3) <sup>4</sup> was used to prepare the reference genome (`bismark_genome_preparation`) from the Sus scrofa Sscrofa11.1 genome<sup>4, 5</sup>. Bismark (0.22.3) <sup>4</sup> was also used for aligning reads to the reference (`bismark`) and extracting the methylation information (`bismark_methylation_extractor`). Alignments were ran with the `pbat` flag to properly align reads to the complementary-to-original strands of the genome and methylation extraction was ran with the `--comprehensive` and `-merge_non_CpG` flags to extract methylation in CpG context, `--bedGraph` to get an output in bedgraph format, and `--ignore_3prime 6` to ignore methylation biases on the 3' ends of the reads. Loci were filtered to include only loci that had at least one count (methylated or unmethylated) in every sample. Methylated and unmethylated counts were summed across Ensembl gene promoter regions (within 2kb of a transcription start site) using biomaRt(2.46.3)<sup>7, 8 </sup> and GRanges(1.42.0)<sup>9</sup> (`distanceToNearest`). The edgeR<sup>10</sup> package was used to find differentially methylated promoters. Briefly, the glmFit function was used to fit a negative binomial generalized log-linear model. The experimental design matrix was constructed using modelMatrixMeth with a factorial experimental design (~0 + group), where group was a factor variable with levels comprised of each combination of treatment/diet/tissue. We dropped the intercept from our model to parameterize it as a means model. Subsequently, the glmLRT function was used to find differentially methylated promoters for comparisons of interest, which were made by constructing contrast vectors. Promoters were considered differentially methylated if the nominal p-value was < 0.01. 

1.) Andrews (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc

2.) Krueger F. (2012). A wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files, with some extra functionality for MspI-digested RRBS-type (Reduced Representation Bisufite-Seq) libraries. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

3.) Bolger et al. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. DOI: 10.1093/bioinformatics/btu170 

4.) Krueger (2010). Bismark: A tool to map bisulfite converted sequence reads and determine cytosine methylation states. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/bismark/

5.) Li et al. (2016). Comprehensive variation discovery and recovery of missing sequence in the pig genome using multiple de novo assemblies. DOI: 10.1101/gr.207456.116 

6.) Warr et al. (2020). An improved pig reference genome sequence to enable pig genetics and genomics research. DOI: 10.1093/gigascience/giaa051

7.) Durinck et al. (2009). Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. DOI: 10.1038/nprot.2009.97 

8.) Durinck et al. (2005). BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis. DOI: 10.1093/bioinformatics/bti525 

9.) Lawrence et al. (2013). Software for Computing and Annotating Genomic Ranges. DOI: 10.1371/journal.pcbi.1003118

10.) Chen et al. (2018).  Differential methylation analysis of reduced representation bisulfite sequencing experiments using edgeR. DOI: 10.12688/f1000research.13196.2 

# Read QC



## Part 7: Methylation Extracor 

Bismark methylation extractor was used to extract methylation information using the following parameters:
```
bismark_methylation_extractor --bedGraph --comprehensive --ignore_3prime 6 -s --merge_non_CpG --report --output /gpfs/data/cbc/fedulov_alexey/porcine_rrbs/alignments/ --gzip --multicore 8 --genome_folder /gpfs/data/shared/databases/refchef_refs/S_scrofa/primary/bismark_index/ $sample  
```
This produced m-bias plots, illustrating the methylation proportion across each possible position in the read (see Part 8 below). 

## Part 8: M-bias plots 
 **Library** | **M-bias plot** | 
|:------:|:-----------:|
| Rh1_R1_001 | ![](png_mbias/rh1_Bismark M-bias Read 1.png) | 
| Rh2_R1_001 | ![](png_mbias/rh2_Bismark M-bias Read 1.png) | 
| Rh3_R1_001 | ![](png_mbias/rh3_Bismark M-bias Read 1.png) | 
| Rh4_R1_001 | ![](png_mbias/rh4_Bismark M-bias Read 1.png) | 
| Rh5_R1_001 | ![](png_mbias/rh5_Bismark M-bias Read 1.png) | 
| Rh6_R1_001 | ![](png_mbias/rh6_Bismark M-bias Read 1.png) | 
| Rh7_R1_001 | ![](png_mbias/rh7_Bismark M-bias Read 1.png) | 
| Rh8_R1_001 | ![](png_mbias/rh8_Bismark M-bias Read 1.png) | 
| Rh10_R1_001 | ![](png_mbias/rh10_Bismark M-bias Read 1.png) | 
| Rh11_R1_001 | ![](png_mbias/rh11_Bismark M-bias Read 1.png) | 
| Rh12_R1_001 | ![](png_mbias/rh12_Bismark M-bias Read 1.png) | 
| Rh14_R1_001 | ![](png_mbias/rh14_Bismark M-bias Read 1.png) | 
| Rh15_R1_001 | ![](png_mbias/rh15_Bismark M-bias Read 1.png) | 
| Rh16_R1_001 | ![](png_mbias/rh16_Bismark M-bias Read 1.png) | 
| Rh18_R1_001 | ![](png_mbias/rh18_Bismark M-bias Read 1.png) | 
| Rh19_R1_001 | ![](png_mbias/rh19_Bismark M-bias Read 1.png) | 
| Rh21_R1_001 | ![](png_mbias/rh21_Bismark M-bias Read 1.png) | 
| Rh22_R1_001 | ![](png_mbias/rh22_Bismark M-bias Read 1.png) | 


## Part 9: Filtering and GLM (Statistical Analysis)
Genes were filtered such that only CpGs with at least 5 counts (methylated and unmethylated) in 3 samples were included for downstream analysis. 
Additionally, CpGs that were never methylated or always methylated were filtered out, as these provide no information about differential methylation. 
The resulting sample size after filtering was 84,999. 

The glmFit function was used to fit a negative binomial generalized log-linear model. The experimental design matrix was constructed using modelMatrixMeth with a factorial experimental design (~0 + group), where group was a factor variable with levels comprised of each combination of treatment/diet/tissue. We dropped the intercept from our model to parameterize it as a means model. Subsequently, the glmLRT function was used to find differentially methylated loci for comparisons of interest, which were made by constructing contrast vectors. Individual CpG sites were considered differentially methylated if the nominal p-value was < 0.01. Results were filtered based on this nominal p-value threshold.  

Results from the statistical analyses have been stored in excel files (and sent to you via email), with files being named by research question. Note that the excel files entitled `both_HSM_and SMV.xlsx`, `SMV_effect_only.xlsx`, and `HSMV_effect_only.xlsx` are the overlap in results lists from the HSMV normal vs HSMV ischemic and SMV normal vs. SMV ischemic comparisons.`both_HSM_and SMV.xlsx` includes loci for which there was an effect for tissue type in both high fat and normal diets; `SMV_effect_only.xlsx` includes loci for which there was an effect for tissue type in normal diet only (and not in high fat); and conversely, `HSMV_effect_only.xlsx` includes loci for which there was an effect of tissue type in high fat diet only (and not in normal diet).

The results files contain the following columns: 

- "`loci` column is the location of a given genomic location. The first number or letter is the chromosome, followed by a `-` symbol and a second number indicates the position or location on the chromosome."

- "`logFC` indicates the log2-fold change of expression between the two conditions tested. For `SMVisch_vs_SMVnormal_results_filtered` results, positive values indicate higher methylation rates in SMVisch relative to SMVnormal, negative values indicate lower methylation rates in SMVisch vs SMVnormal. For `HSMVisch_vs_HSMVnormal_results_filtered` results, positive values indicate higher methylation rates in HSMVisch relative to HSMVnormal, negative values indicate lower methylation rates in HSMVisch relative to HSMVnormal. For `diffisch_vs_diffnormal_results_filtered` results, positive values indicate that the difference due to normal vs ischemic (i.e., effect of ischemia) is greater in HSMV than SMV, whereas negative values indicate that the difference due to normal vs ischemic (i.e., effect of ischemia) is greater in SMV than in HSMV. For `MVM_vs_SMVischemic_results_filtered` results, positive values indicate higher methylation rates in MVM relative to SMVischemic, negative values indicate lower methylation rates in MVM relative to SMVischemic." Lastly, for `HVM_vs_HSMVischemic_results_filtered` results, positive values indicate higher methylation rates in HVM relative to HSMVischemic, negative values indicate lower methylation rates in HVM relative to HSMVischemic."

- "`logCPM` is the average log2-counts per million, the average taken over all libraries used in the experiment."

- "`LR` is the likelihood ratio statistics (larger LR, smaller p-value)"

- "`PValue` probability of obtaining results at least as extreme as the results actually observed, under the assumption that the null hypothesis is correct"

- "`BH` is Benjamini & Hochberg corrected p-values"

- "`bonferroni` Bonferroni corrected p-values"

- "The columns with the format `#-Me` or `#-Un` indicate the methylated and unmethylated counts at a given locus for a particular sample."

- "`Chr` is the chromosome of the locus"

- "`Locus` is the position of the locus on the chromosome indicated in the `Chr` column."

- "The columns with the format `methylated_proportion_#` indicate the proportion of methylated counts at a given loci for a particular sample, calculated as `methylated counts / (unmethylated counts + methylated counts)` for each sample and loci."

- "`TSS_start` indicates the TSS location nearest to the methylation loci shown."

- "`uniprot_gn_symbol` is the gene symbol (if available) associated with each TSS"

- "`ensembl_gene_id` is the Ensembl gene ID associated with each TSS"

- "`gene_start` and `gene_end` are the start and end genomic positions for each gene."




## Part 2: Read quality plots 
 **Library** | **Raw reads**  | **Trimmed reads** |
|:------:|:-----------:|:----------:|
| Rh1_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh1_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh1_trim_png.png) |
| Rh2_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh2_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh2_trim_png.png) |
| Rh3_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh3_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh3_trim_png.png) |
| Rh4_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh4_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh4_trim_png.png) |
| Rh5_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh5_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh5_trim_png.png) |
| Rh6_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh6_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh6_trim_png.png) |
| Rh7_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh7_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh7_trim_png.png) |
| Rh8_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh8_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh8_trim_png.png) |
| Rh10_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh10_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh10_trim_png.png) |
| Rh11_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh11_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh11_trim_png.png) |
| Rh12_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh12_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh12_trim_png.png) |
| Rh14_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh14_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh14_trim_png.png) |
| Rh15_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh15_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh15_trim_png.png) |
| Rh16_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh16_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh16_trim_png.png) |
| Rh18_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh18_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh18_trim_png.png) |
| Rh19_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh19_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh19_trim_png.png) |
| Rh21_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh21_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh21_trim_png.png) |
| Rh22_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/Rh22_png.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/Rh22_trim_png.png) |

## Part 3: Per base sequence content
**Library** | **Raw reads**  | **Trimmed reads** |
|:------:|:-----------:|:----------:|
| Rh1_R1_001 | ![](base/Rh1_base.png) | ![](base_trim/Rh1_base_trim.png) |
| Rh2_R1_001 | ![](base/Rh2_base.png) | ![](base_trim/Rh2_base_trim.png) |
| Rh3_R1_001 | ![](base/Rh3_base.png) | ![](base_trim/Rh3_base_trim.png) |
| Rh4_R1_001 | ![](base/Rh4_base.png) | ![](base_trim/Rh4_base_trim.png) |
| Rh5_R1_001 | ![](base/Rh5_base.png) | ![](base_trim/Rh5_base_trim.png) |
| Rh6_R1_001 | ![](base/Rh6_base.png) | ![](base_trim/Rh6_base_trim.png) |
| Rh7_R1_001 | ![](base/Rh7_base.png) | ![](base_trim/Rh7_base_trim.png) |
| Rh8_R1_001 | ![](base/Rh8_base.png) | ![](base_trim/Rh8_base_trim.png) |
| Rh10_R1_001 | ![](base/Rh10_base.png) | ![](base_trim/Rh10_base_trim.png) |
| Rh11_R1_001 | ![](base/Rh11_base.png) | ![](base_trim/Rh11_base_trim.png) |
| Rh12_R1_001 | ![](base/Rh12_base.png) | ![](base_trim/Rh12_base_trim.png) |
| Rh14_R1_001 | ![](base/Rh14_base.png) | ![](base_trim/Rh14_base_trim.png) |
| Rh15_R1_001 | ![](base/Rh15_base.png) | ![](base_trim/Rh15_base_trim.png) |
| Rh16_R1_001 | ![](base/Rh16_base.png) | ![](base_trim/Rh16_base_trim.png) |
| Rh18_R1_001 | ![](base/Rh18_base.png) | ![](base_trim/Rh18_base_trim.png) |
| Rh19_R1_001 | ![](base/Rh19_base.png) | ![](base_trim/Rh19_base_trim.png) |
| Rh21_R1_001 | ![](base/Rh21_base.png) | ![](base_trim/Rh21_base_trim.png) |
| Rh22_R1_001 | ![](base/Rh22_base.png) | ![](base_trim/Rh22_base_trim.png) |

## Part 4: GC content plots 
 **Library** | **Raw reads**  | **Trimmed reads** |
|:------:|:-----------:|:----------:|
| Rh1_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh1_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh1_trim_gc.png) |
| Rh2_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh2_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh2_trim_gc.png) |
| Rh3_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh3_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh3_trim_gc.png) |
| Rh4_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh4_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh4_trim_gc.png) |
| Rh5_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh5_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh5_trim_gc.png) |
| Rh6_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh6_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh6_trim_gc.png) |
| Rh7_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh7_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh7_trim_gc.png) |
| Rh8_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh8_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh8_trim_gc.png) |
| Rh10_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh10_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh10_trim_gc.png) |
| Rh11_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh11_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh11_trim_gc.png) |
| Rh12_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh12_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh12_trim_gc.png) |
| Rh14_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh14_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh14_trim_gc.png) |
| Rh15_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh15_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh15_trim_gc.png) |
| Rh16_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh16_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh16_trim_gc.png) |
| Rh18_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh18_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh18_trim_gc.png) |
| Rh19_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh19_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh19_trim_gc.png) |
| Rh21_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh21_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh21_trim_gc.png) |
| Rh22_R1_001 | ![](fastqc_reports_Alexey/initial_qc/png/rh22_gc.png) | ![](fastqc_reports_Alexey/retrim_qc/png_trim/rh22_trim_gc.png) |



## Part 6: Alignment summary statistics
| **LIBRARY** | **NUMBER_READS** | **MAPPED_READS** | **UNMAPPED_READS** | **DUPLICATED_READS** | **DUPLICATION_RATE** |  
|:---:|:---:|:---:|:---:|:---:|:---:|
| Rh1_R1_001 | 10,845,489 | 10,845,489 / 100% | 0 / 0% | 1,437,581 / 13.26% | 12.6% | 
| Rh2_R1_001 | 12,892,543 | 12,892,543 / 100% | 0 / 0% | 1,496,670 / 11.61% | 10.93% | 
| Rh3_R1_001 | 8,715,481 | 8,715,481 / 100% | 0 / 0% | 1,029,908 / 11.82% | 11.02% | 
| Rh4_R1_001 | 9,968,480 | 9,968,480 / 100% | 0 / 0% | 1,209,925 / 12.14% | 11.51% | 
| Rh5_R1_001 | 13,925,551 | 13,925,551 / 100% | 0 / 0% | 2,056,379 / 14.77% | 13.91% | 
| Rh6_R1_001 | 13,362,404 | 13,362,404 / 100% | 0 / 0% | 1,750,619 / 13.1% | 12.31% | 
| Rh7_R1_001 | 12,851,559 | 12,851,559 / 100% | 0 / 0% | 1,713,430 / 13.33% | 12.48% | 
| Rh8_R1_001 | 12,872,887 | 12,872,887 / 100% | 0 / 0% | 1,769,849 / 13.75% | 12.93% | 
| Rh10_R1_001 | 15,060,158 | 15,060,158 / 100% | 0 / 0% | 2,107,684 / 14% | 13.11% | 
| Rh11_R1_001 | 13,432,863 | 13,432,863 / 100% | 0 / 0% | 1,930,952 / 14.37% | 13.51% | 
| Rh12_R1_001 | 14,880,989 | 14,880,989 / 100% | 0 / 0% | 2,381,186 / 16% | 14.75% | 
| Rh14_R1_001 | 13,023,116 | 13,023,116 / 100% | 0 / 0% | 1,938,398 / 14.88% | 13.62% | 
| Rh15_R1_001 | 13,292,288 | 13,292,288 / 100% | 0 / 0% | 1,992,503 / 14.99% | 13.72% | 
| Rh16_R1_001 | 17,058,927 | 17,058,927 / 100% | 0 / 0% | 2,756,384 / 16.16% | 15.03% | 
| Rh18_R1_001 | 12,196,440 | 12,196,440 / 100% | 0 / 0% | 1,814,764 / 14.88% | 13.42% | 
| Rh19_R1_001 | 13,248,129 | 3,248,129 / 100% | 0 / 0% | 1,901,092 / 14.35% | 13.16% | 
| Rh21_R1_001 | 16,030,293 | 16,030,293 / 100% | 0 / 0% | 2,624,856 / 16.37% | 14.84% | 
| Rh22_R1_001 | 13,927,765 | 13,927,765 / 100% | 0 / 0% | 1,961,472 / 14.08% | 13.05% | 

| **LIBRARY** | **COVERAGE_MEAN** | **COVERAGE_SD** | **MEAN_MAP_QUALITY** | 
|:---:|:---:|:---:|:---:|
| Rh1_R1_001 | 0.1844 | 0.7664 | 28.42 |
| Rh2_R1_001 | 0.2191 | 0.8243 | 29.32 | 
| Rh3_R1_001 | 0.1482 | 0.7405 | 27.77 |
| Rh4_R1_001 | 0.1694 | 0.7166 | 28.62 | 
| Rh5_R1_001 | 0.2368 | 0.943 | 28.68 | 
| Rh6_R1_001 | 0.2271 | 0.9149 | 29.01 | 
| Rh7_R1_001 | 0.2185 | 0.8851 | 28.47 | 
| Rh8_R1_001 | 0.2188 | 0.9475 | 29.15 | 
| Rh10_R1_001 | 0.2561 | 1.0472 | 28.59 | 
| Rh11_R1_001 | 0.2284 | 0.9772 | 28.47 | 
| Rh12_R1_001 | 0.253 | 1.0263 | 29.06 | 
| Rh14_R1_001 | 0.2214 | 0.9804 | 29.35 | 
| Rh15_R1_001 | 0.226 | 0.9128 | 28.77 | 
| Rh16_R1_001 | 0.29 | 1.2004 | 29.18 | 
| Rh18_R1_001 | 0.2075 | 1.0798 | 28.59 | 
| Rh19_R1_001 | 0.2252 | 0.9003 | 29.01 | 
| Rh21_R1_001 | 0.2726 | 1.2379 | 29.02 | 
| Rh22_R1_001 | 0.2368 | 1.1062 | 29.13 | 