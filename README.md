# Models of trinucleotide retention in reverse phased liquid chromatography
This repository is a sumplement information of article: "Thou shalt know thy neighbour – A sequence dependent separation of trinucleotides by ion-interaction reversed-phase liquid chromatography" by Mikulášek et al.. 
All scripts used to compute pushlished results are present. 
The organisation is as follows:

./data/ - folder of all input data: simulation series of SASA and egb, retention times

./scripts/ - all computing R scripts

./execute_all.R - calls all computing scripts on apripriate data and saves results to tables in this folder

./README.md - this readme.

R packages DAAG and MASS are needed to compute the output. If they are present, the following command should execute all scripts:

Rscript ./execute_all.R
