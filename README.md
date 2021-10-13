# Fedulov Porcine RRBS

This repository tracks the analysis pipeline for the RRBS analysis of porcine data for the Fedulov lab. All analysis steps are recorded in the `Alexey_RRBS_Porcine_data_analysis_pipeline.md` markdown file. The scripts folder contains all of the scripts used to run various aspects of the analysis. The portions of `Alexey_RRBS_Porcine_data_analysis_pipeline.md` that use bash were run on Oscar, while the R sections were run on a local machine.

To build the docker container:

```{bash}
docker build -t rrbs:latest .
```

To tag and push to dockerhub:

```{bash}
docker tag rrbs compbiocore/rrbs:latest
docker push compbiocore/rrbs:latest
```

To run the container, make sure you are in the folder you want to mount, then run::

```{bash
docker run --rm -p 8787:8787 -e USER=rstudio -e PASSWORD=yourpassword --volume ${PWD}:/home/rstudio rrbs:latest
```

Then navigate to localhost:8787 in firefox or chrome.


To run as Singularity on Oscar:

```{bash}
cd /gpfs/data/cbc/fedulov_alexey/porcine_rrbs/singularity
mkdir -p run var-lib-rstudio-server
printf 'provider=sqlite\ndirectory=/var/lib/rstudio-server\n' > database.conf
singularity pull docker://compbiocore/rrbs:latest
```

Then log into oscar over VNC client, then open terminal and run:

```{bash}
cd /gpfs/data/cbc/fedulov_alexey/porcine_rrbs/singularity
export SINGULARITY_BINDPATH="/gpfs/data/cbc"
singularity exec --bind run:/run,var-lib-rstudio-server:/var/lib/rstudio-server,database.conf:/etc/rstudio/database.conf,gpfs/data/cbc/fedulov_alexey/porcine_rrbs:/home/rstudio/porcine_rrbs compbiocore-rrbs.sif rserver --www-address=127.0.0.1
```

Then open another terminal window, load firefox module, and navigate to localhost:8787
