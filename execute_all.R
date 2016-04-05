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
print("Sequence dependent retention of trinucleotides were computed, summary was saved to ./output/1_BSmodel.txt and ./output/1_ASmodel.txt recpectively")

# ---- permutation test of AS-model, warning, takes few minutes ----
source('./scripts/s1_permutation_test.R')
print("Distribution of p values and r squared of randomized BS a AM model were saved to ./output/s1_...")

# ---- sasa NNN-model ----
source('./scripts/2_egb_sasa.R')
print('model_NNN using egb and sasa (equation 10) was constructed, summary was saved to ./output/2_NNN_model_summary.txt')

# ---- detection deviant nucleoties in NNN-model ----
source('./scripts/3_deviant_nucleotides.R')
print('tests for the most deviant nucleotides (equation 11, table 3) , table and summary was saved to ./output/3_deviant_nucleotides.csv and ./output/3_NGA_model_summary.txt respectively')

# ---- sasa standard deviatino estimation ----
source('./scripts/s2_sasa_sd_estimation.R')
print('Aproxximation of standard deviation of sasa model, output was saved to ./output/s2_sasa_sd.csv')



