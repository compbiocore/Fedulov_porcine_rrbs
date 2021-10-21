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
singularity exec --bind run:/run,var-lib-rstudio-server:/var/lib/rstudio-server,database.conf:/etc/rstudio/database.conf,/gpfs/data/cbc/fedulov_alexey/Fedulov_porcine_rrbs:/home/rstudio/Fedulov_porcine_rrbs compbiocore-rrbs.sif rserver --www-address=127.0.0.1
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
6. and “HSMV ischemic” = myocardial injection of vehicle control in a high-fat diet – ischemic heart tissue.

Goal: we want to determine the loci differentially methylated in all pairwise comparisons between the groups to answer these experimental questions:

1). What was the effect of ischemia vs. normal heart tissue on a normal diet in absence of treatments? (SMV normal vs. SMV ischemic).

2). What was the effect of ischemia vs. normal heart tissue on a high-fat diet in absence of treatments (HSMV normal vs. HSMV ischemic) and whether it was different from that effect in the normal diet? (this resultant list vs. list from #1) 

3). What was the effect of microvesicle treatment (vs. ischemia) in normal diet and in high-fat diet? (MVM vs. SMV ischemic; and HVM vs. HSMV ischemic)

Notes: The differentially methylated loci detection should be done at reasonable stringency to a) characterize the extent of changes by each experimental factor, and b) to inform downstream pathway analysis.
The differentially methylated loci need to be mapped to the nearest transcript (with distance to TSS) and presented in the form of data tables convenient for heatmaps (i.e. with methylation percent values) and for downstream pathway analysis (i.e. with gene names and IDs).
