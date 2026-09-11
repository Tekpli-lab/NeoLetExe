This repository accompanies the study "Longitudinal Single-Cell Profiling of Breast Tumors Reveals Immune and Epigenetic Mechanisms of Sensitivity and Resistance to Aromatase Inhibitors", published in XXXXXX. Its purpose is to reproduce key figures and analyses from the paper.

**Aim of the Study and Main Results**

We analyzed scRNA-seq, scTCR-seq, and scBCR-seq data from breast tumor biopsy specimens collected at different time points during neoadjuvant therapy with aromatase inhibitors. This study investigates tumor evolution under treatment pressure and deciphers the mechanisms of sensitivity and resistance to hormone therapy.

**Repository Structure**

- **pre-processing**: Contains scripts to obtain all the key input files needed to produce the figures and describes the logic of the data analysis.
- **Figures** and **Supplementary_Figures**: Contains code to plot all figures and supplementary figures in the paper.
- **ShinyApp**: Provides a user-friendly interface for interacting with the data, allowing for interactive analyses of T/NK cells, myeloid cells, and all cell UMAPs without any coding required.

**Data Access**

T/NK cells, myeloid cells, and all cell UMAPs are available for exploration a ShinyApp. 
Processed data are accessible at:  https://ous-research.no/home/tekpli/Research-projects
Raw CellRanger scRNA-seq, scTCR-seq, and scBCR-seq data have been deposited in the European Genome-Phenome Archive (EGA), hosted by the EBI and the CRG, under the accession code EGAS00000002766.


[![DOI](https://zenodo.org/badge/1069113536.svg)](https://doi.org/10.5281/zenodo.22704982)
