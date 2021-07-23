# Fedulov Porcine RRBS

This repository tracks the analysis pipeline for the RRBS analysis of porcine data for the Fedulov lab. All analysis steps are recorded in the Alexey_RRBS_Porcine_data_analysis_pipeline.md markdown file. The fedulov_rrbs folder contains a scripts subfolder, which contains all of the scripts used to run various aspects of the analysis. The portions of analysis.md that use bash were run on Oscar, while the R sections were run on a local machine.

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
