This directory contains the code to run a snakemake pipeline for doing the analysis of cis-regulatory elements on SCAFE counts. 

In order to run, you will need to have snakemake and singularity installed. Dependencies are handeled through containers downloaded automatically from dockerhub.
In the config.yaml you will get an overview of what input files are needed to run the pipeline. 

When the input data is in order, run the snakemake pipeline from the current directory with this command, N being the number of cores you wish to allocate:
snakemake --use-singularity -c N
