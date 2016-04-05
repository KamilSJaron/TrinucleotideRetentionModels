# Models of trinucleotide retention in reverse phased liquid chromatography
This repository is a sumplement information to paper XXXXXX, showing exact script used to acess pushlished results. 
The organisation is as follows:

./data/ - folder of all input data: SASA monte carlo simulations and retention times (will be added in the date of publicaiton)

./scripts/ - all computing R scripts

./execute_all.R - calls all computing scripts on apripriate data and saves results to tables in this folder

./README.md - this readme.

R packages DAAG and MASS are needed to compute the output. If they are present, the following command should execute all scripts:

Rsrcipt ./execute_all.R
