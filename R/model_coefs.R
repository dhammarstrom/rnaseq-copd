##-------------------------------------
## model_coef.R
##
## Title: Extract model coefficients from mixed-model output
## Purpose: A single file to get a data frame of model coefficients from mixed models
## Author:
##
##
##
##-------------------------------------
## Notes:
# This script stores a data frame with model coefficients in the environment. 
# Filtering based on uniformity is done here.
#
#
#
#
#
#
## ------------------------------------


source("./R/lib_fun.R")


# All models and coefs 
seq_data <- readRDS(file = "./data/derivedData/mixed_model/mixed_full_copd.RDS")


# Load all coefs and store in environment  

all_coefs <- seq_data %>%
  group_by(coef, model) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  ungroup() %>%
  filter(uniform.p > 0.05) %>%
  dplyr::select(gene, model, interaction, coef, estimate, se, z.val, p.val, p.adj) 
  




