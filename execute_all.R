# clear all variables
rm(list = ls())

# ---- Packages ----

library(MASS)
library(DAAG)

# ---- Load data ----
source('./scripts/0_load_data.R')
print("data were loaded")

# ---- Linear regression model of sequence dependent retention of trinucleotides ----
source('./scripts/1_sequence_dependent_models.R') 
print("equence dependent retention of trinucleotides were computed, summary was saved to ./output/1_BSmodel.csv and ./output/1_ASmodel.csv recpectively")

# ---- permutation test of AS-model ----
source('./scripts/s1_permutation_test.R')
print("Distribution of p values and r squared of randomized BS a AM model were saved to ./output/s1_...")

# ---- sasa NNN-model ----
source('./scripts/2_egb_sasa.R')
print('Model using ???? and sasa (equation 10) was constructed, summary was saved to ./output/2_SASAmodel.csv')

# ---- detection deviant nucleoties in NNN-model ----
source('./scripts/3_deviant_nucleotides.R')
print('XXXXXXXXX (equation 11, table 3) , summary was saved to ./output/3_deviant_nk.csv')

# ---- sasa standard deviatino estimation ----
source('./scripts/s2_sasa_sd_estimation.R')
print('XXXXXXXXX , summary was saved to ./output/s2_sasa_sd.csv')