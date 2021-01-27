##-------------------------------------
## edgeR_presupp.R
##
## Title:  EdgeR DE analysis
## Purpose: 
## Author:
##
##
##
##-------------------------------------
## Notes: Expression differences between COPD and control at baseline
# Using edgeR as the design at baseline does not include repeated measures. 
#
#
#
#
#
#
#
## ------------------------------------

source("./R/lib_fun.R")




# Read dge list
dge <- readRDS(file = "./data/derivedData/dge_lists/dge_sum.RDS")

# Muscle weight data
dat_yield <- read_excel("./data/subject_sample.xlsx", na = "NA") %>%
  inner_join(read_excel("./data/subject_totalrna.xlsx", na = "NA") %>%
               mutate(RM_leg = as.character(RM_leg))) %>%
  mutate(time = if_else(time == "3W", "W3", time)) %>%
  dplyr::select(ex.nr, yield, kolsfrisk) %>%
  mutate(ex.nr = as.character(ex.nr), 
         copd = if_else(kolsfrisk == 0, "control", "copd"), 
         copd = factor(copd, levels = c("control", "copd")))


## Read sex
subject_sex <- read_excel("./data/subject_sex.xlsx") %>%
  mutate(sex = if_else(sex == 0, "female", "male"), 
         subject = as.character(subject)) 




# Subset the dge list to only include PreSupp
dge_presupp <- dge[,rownames(dge$samples) %in% rownames(dge$samples[dge$samples$time == "PreSupp",]) ]

# Add information (sex and copd status) to sample data frame
dge_presupp$samples <- dge_presupp$samples %>%
  inner_join(subject_sex) %>%
  inner_join(dat_yield) %>%
  print()

# Fix the deisgn matrix
dge_presupp$design <- model.matrix(~ sex + copd, data = dge_presupp$samples)


# Estimate dispersion 
dge_presupp <- estimateDisp(dge_presupp, design = dge_presupp$design, 
                            robust = TRUE)

# Fit models
fitQL <- glmQLFit(dge_presupp, design = dge_presupp$design)

# Perform test
f_test <- glmQLFTest(fitQL, coef = "copdcopd", poisson.bound = TRUE)

# Retrieve estimates and calculate CI 
results_edgeR <- topTags(f_test, n = Inf) %>%
  data.frame() %>%
  mutate(gene = rownames(.), 
         z = -0.862 + sqrt((0.743-2.404 * log(PValue))), 
         se = logFC/z, 
         ciu = logFC + se * 1.96, 
         cil = logFC - se * 1.96) %>%
  print()

# Save results
saveRDS(results_edgeR, file = "./data/derivedData/edgeR/presupp_results.RDS")
