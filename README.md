# Analysis from the "Simulation of 69 microbial communities indicates sequencing depth and false positives are major drivers of bias in Prokaryotic metagenome-assembled genome recovery." study.

## Introduction
This repository contains the R code and data for the research conducted in [Title of Your Paper.............]. 

In this study, we tackle MAG recovery knowledge gaps by generating mock communities with varying parameters of sequencing depth, taxonomic distribution relatedness and species abundance profiles. We use these datasets to evaluate three pipelines for binning, each using a different approach to bin the sequences, and assess which factors or combinations of factors drive Prokaryotic MAG recovery and if these are pipeline-dependentpipeline dependent. Additionally, we evaluate the effect of multiple or linear chromosomes in MAG recovery.

For a more detailed explanation of the study, please refer to the original publication
## Citation

Nunes da Rocha, U., Coelho Kasmanas, J., Toscan, R., Sanches, D. S., Magnusdottir, S., & Saraiva, J. (2023). Simulation of 69 microbial communities indicates sequencing depth and false positives are major drivers of bias in Prokaryotic metagenome-assembled genome recovery. bioRxiv, 2023-05.

## Repository Contents

### Scripts

The running of this repository is meant to happen sequentially in numerical folder order. Inside each folder, the user should find an R script also properly ordered by number, and the required input data. The input data is either generated in a previous step of the analysis pipeline (this repository) or was generated in the original-publication described manner. 

### Running the Scripts
To reproduce the results or to conduct further analysis, the user should have a working installation of R version 4.2.1 (although it may also work on previous versions) and the required R packages loaded at the initial run of the script.
The required R packages are:

- tidyr v1.2.1
- ggplot2 v3.4.0
- dplyr v1.0.10
- ggpubr v0.5.0
- rstatix v0.7.1


---


**Acknowledgments**: 

[Any acknowledgments to contributors, institutions, funding sources, etc.]

