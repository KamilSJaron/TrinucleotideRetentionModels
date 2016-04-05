# clear all variables
rm(list = ls())

# ---- Packages ----

library(MASS)
library(DAAG)

# ---- Load data ----
source('./scripts/0_load_data.R')

# ---- sequential models ----

# ---- permutation test of AS-model ----
source('./scripts/s1_permutation_test.R')

