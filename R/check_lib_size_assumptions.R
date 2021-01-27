##-------------------------------------
## check_lib_size_assumptions.R
##
## Title:
## Purpose: Check assumptions of equal tissue/tot-RNA 
## Author:
##
##
##
##-------------------------------------
## Notes: To check if any group differences existed for library size, tissue used in cDNA prep
# or library size per tissue, these variables were tested over the course of the study. 
# There were no indications of (strong) systematic differences between groups. 
# Lib-size as normalization is the preferred model as tissue offset would not give more information 
# on the group level.
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
         kolsfrisk = if_else(kolsfrisk == 0, "control", "copd"), 
         kolsfrisk = factor(kolsfrisk, levels = c("control", "copd")))


## Read sex
subject_sex <- read_excel("./data/subject_sex.xlsx") %>%
  mutate(sex = if_else(sex == 0, "female", "male"), 
         subject = as.character(subject)) 




dat <- dat_yield %>%
  inner_join(dge$samples) %>%
  filter(time != "PreSupp") %>%
  mutate(time = factor(time, levels = c("PreExc", "W3", "PostExc")), 
         tissue = 1000/yield) %>%
  
  print()



# RNA-yield per tissue COPD vs. control

m <- lme(log(yield) ~ time * kolsfrisk, random = list(subject = ~ 1, 
                                                 leg = ~1), 
         data = dat)

plot(m)
# OK resid plot 

# Extract marginal means
emmeans(m, specs = ~time | kolsfrisk) %>%
  data.frame() %>%
  ggplot(aes(time, exp(emmean), color = kolsfrisk, group = kolsfrisk)) + geom_line() + 
  geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), 
                width = 0.2)

# Similar patterns in both groups
# Confirmed that interaction effect on group level was negliable
anova(m)
summary(m)


# Library size tissue, and lib size per tissue between groups ##############

lib.m <- lme(log(lib.size * norm.factors) ~ time * kolsfrisk, random = list(subject = ~ 1, 
                                                               leg = ~1), 
                  data = dat)

lib.tissue.m <- lme(log((lib.size * norm.factors)/tissue) ~ time * kolsfrisk, random = list(subject = ~ 1, 
                                                                            leg = ~1), 
             data = dat)


tissue.m <- lme(log(tissue) ~ time * kolsfrisk, random = list(subject = ~ 1, 
                                                                                            leg = ~1), 
                    data = dat)





plot(lib.tissue.m)
plot(lib.m)
plot(tissue.m)
# OK resid plot 

# Extract marginal means
emmeans(lib.m, specs = ~time | kolsfrisk) %>%
  data.frame() %>%
  ggplot(aes(time, exp(emmean), color = kolsfrisk, group = kolsfrisk)) + geom_line() + 
  geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), 
                width = 0.2)


# Extract marginal means
emmeans(lib.tissue.m, specs = ~time | kolsfrisk) %>%
  data.frame() %>%
  ggplot(aes(time, exp(emmean), color = kolsfrisk, group = kolsfrisk)) + geom_line() + 
  geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), 
                width = 0.2)


emmeans(tissue.m, specs = ~time | kolsfrisk) %>%
  data.frame() %>%
  ggplot(aes(time, exp(emmean), color = kolsfrisk, group = kolsfrisk)) + geom_line() + 
  geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), 
                width = 0.2)



anova(tissue.m)
anova(lib.tissue.m)
anova(lib.m)



summary(tissue.m)
summary(lib.tissue.m)
summary(lib.m)





