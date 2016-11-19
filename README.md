# Models of trinucleotide retention in reverse phased liquid chromatography

This repository contains all the `R` scripts used for calculations in article: "Sequence-dependent separation of trinucleotides by ion-interaction reversed-phase liquid chromatography—A structure-retention study assisted by soft-modelling and molecular dynamics" by Mikulášek et al., DOI: [10.1016/j.chroma.2016.09.060](http://www.sciencedirect.com/science/article/pii/S0021967316312766). 

The organisation is as follows:

[./data/](data) - folder of all input data: simulation series of SASA and egb, retention times

[./scripts/](scripts) - all computing `R` scripts

[./execute_all.R](execute_all.R) - calls all computing scripts on apripriate data and saves results to tables in this folder

./README.md - this readme.

`R` packages `DAAG` and `MASS` are required. If they are installed, the following command should execute all scripts:

```bash
Rscript ./execute_all.R
```
