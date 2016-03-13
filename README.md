# Models of trinucleotide retention in reverse phased liquid chromatography
This repository is a sumplement information to paper XXXXXX, showing exact script used to acess pushlished results. 
The organisation is as follows:

./input_data/ - folder of all input data: SASA monte carlo simulations and retention times

./scripts/ - all computing R scripts

./exectureAll.R - calls all computing scripts on apripriate data and saves results to tables in this folder

./README.md - this readme.

R packages DAAG and MASS are needed to compute the output. If they are present, the following command should execute all scripts:

./Rscipt exectureAll.R
