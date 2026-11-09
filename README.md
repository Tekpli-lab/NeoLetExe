This repository accompanies the study "Longitudinal single-cell profiling of breast tumors highlights immune and epigenetic mechanisms associated with sensitivity to aromatase inhibitors", published in Nature Genetics 4th Nov 2026 (https://doi.org/10.1038/s41588-026-02811-8). Its purpose is to allow to reproduce key figures and analyses.

**Aim of the Study and Main Results**

We analyzed scRNA-seq, scTCR-seq, and scBCR-seq data from breast tumor biopsy specimens collected at different time points during neoadjuvant therapy with aromatase inhibitors. This study investigates tumor evolution under treatment pressure and deciphers the mechanisms of sensitivity and resistance to hormone therapy.

**Repository Structure**

- **pre-processing**: Contains scripts to obtain all the key input files needed to produce the figures and describes the logic of the data analysis.
- **Figures** and **Supplementary_Figures**: Contains code to plot all figures and supplementary figures in the paper.

**Data Access**

T/NK cells, myeloid cells, and all cell UMAPs are available for exploration through a ShinyApp. https://tekpli.shinyapps.io/shinyappmulti_neoletexe_scrnaseq/

Processed data can be downloaded here: https://ous-research.no/home/tekpli/Data-Access
Raw CellRanger scRNA-seq, scTCR-seq, and scBCR-seq data have been deposited in the European Genome-Phenome Archive (EGA), hosted by the EBI and the CRG, under the accession code EGAS00000002766.

Please cite https://doi.org/10.1038/s41588-026-02811-8 when using the code in this repository.
This code has also been deposited on zenodo for archival purposes
[![DOI](https://zenodo.org/badge/1069113536.svg)](https://doi.org/10.5281/zenodo.22704982)
