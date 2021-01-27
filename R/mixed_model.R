### Mixed model approach ##################

light.load <- TRUE # Only load packages needed 
source("./R/lib_fun.R")


# Preliminaries: 
# dge lists are store in ./data/derivedData/dge_lists/dge.RDS
# dge lists are used for gene-wise iterative modeling. 
# The loop is used for one selected list...
# See below for settings on fitting algorithm 
# and diagnostics.


# Notes: 

# RNA-seq data analysis using itearitive fitting of a Poisson- or negative binomial
# generalized mixed model. 
# The main purpose of the analysis is to identify differentially expressed genes between 
# training volume-conditions and time-point regardless of condition. This implementation is inspired by:

# Cui S, Ji T, Li J, Cheng J, Qiu J. 
# What if we ignore the random effects when analyzing RNA-seq data in a multifactor experiment. 
# Stat Appl Genet Mol Biol. 2016;15(2):87–105. doi:10.1515/sagmb-2015-0011

# The normalization factor used in Cui et al:

# "The library scaling factor adjusts for the differential expression 
# caused by both the differential sequencing depths and the differential RNA 
# compositions of different RNA samples. It can be obtained by dividing the 
# effective library size of each library to that of a reference library. 
# Here the effective library size refers to the product of the original 
# library size, the total number of read counts per library, and a normalization 
# factor to adjust for the  RNA  composition  effect.  Throughout  the  paper,  
# we  use  the  trimmed  mean  method  (TMM)  of  Robinson  and Oshlack (2010) 
# to calculate the normalization factor for RNA composition effect, which uses a 
# weighted trimmed  mean  of  the  log  expression  ratios  across  genes  to  
# estimate  the  global  fold  change  of  two  samples  to adjust for the RNA 
# composition effect, assuming that majority of genes are non-differentially expressed."

# edgeR uses effective lib size as an affset in its model (log(eff.lib)), cui suggests 
# using normalization factor for each gene as a fixed factor. These are both fitted below.





# If diagnostics are used, the DHARMa package will test for overdispersion/underdispersion
# and zero inflation 




#### Read data and preliminary exploration ##########



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




##### Approach 2: ################



# Set what dge list to use
selected.dge <- dge

# A list of genes (available in the dge-list)
genes <- rownames(selected.dge)

# set up parallel processing using the foreach package 
cores <- as.integer(round(detectCores() * 0.8, 0))
cl <- makeCluster(cores[1]) # not to overload cpu #### Change this if on server!
registerDoSNOW(cl)

# Progress bar (should work on linux also)
iterations <- length(genes) 
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)



# foreach loop 
results <- foreach(i = 1:length(genes), 
                   .packages = c("glmmTMB", "dplyr", "DHARMa", "car"),  # Packages used in fitting/data wrangling
                   .options.snow = opts) %dopar% {
                           
                           
                           # settings
                           # fitting
                           
                           diagnostics <- TRUE
                           
                           # glmmTMB allows for diagnostics
                           
                         
                           
                           tryCatch(
                                   expr = {
                                           
                                           
                                           # Define which quantile function 
                                           # the function is used for selecting the median sized library
                                           # which in turn is used as reference library as suggested in Cui
                                           
                                           
                                           which.quantile <- function (x, probs = 0.5, na.rm = FALSE){
                                                   if (! na.rm & any (is.na (x)))
                                                           return (rep (NA_integer_, length (probs)))
                                                   
                                                   o <- order (x)
                                                   n <- sum (! is.na (x))
                                                   o <- o [seq_len (n)]
                                                   
                                                   nppm <- n * probs - 0.5
                                                   j <- floor(nppm)
                                                   h <- ifelse((nppm == j) & ((j%%2L) == 0L), 0, 1)
                                                   j <- j + h
                                                   
                                                   j [j == 0] <- 1
                                                   o[j]
                                           }
                                           
                                           
                                           
                                           ## Calculate reference library (this is set as the median library)
                                           reference.lib <- selected.dge$samples[which.quantile(selected.dge$samples$lib.size, na.rm = TRUE),c(2,3)]
                                           
                                           ## Calculate effective library size for the reference library
                                           reference.lib <- reference.lib[1,1] * reference.lib[1,2]
                                           
                                           # Extract data for each sample and calculate normalization factor
                                           
                                      
                                           dat <- selected.dge$samples %>%
                                                   mutate(sample = rownames(.), 
                                                          nf = (lib.size * norm.factors) / reference.lib) 
                                           
                                           ## Extract gene counts and put together in the data frame
                                           dat <- data.frame(counts = selected.dge$counts[genes[i],], 
                                                             sample = colnames(selected.dge$counts)) %>%
                                                   inner_join(dat) %>%
                                                   inner_join(dat_yield) %>%
                                             inner_join(subject_sex) %>%
                                                   filter(time != "PreSupp") %>%
                                                  
                                                   mutate(counts = as.integer(round(counts, 0)), 
                                                          subject = factor(subject), 
                                                          treat = factor(treat, levels = c("placebo", "vitd")), 
                                                          time = factor(time, levels = c("PreExc",
                                                                                         "W3", 
                                                                                         "PostExc")),
                                                          copd = factor(kolsfrisk, levels = c("control", "copd")),
                                                          leg = paste0(leg, subject),
                                                          eff.lib = lib.size * norm.factors, 
                                                          tissue = log(1000/yield))
                                           
                                           
                                           
                                           # A list of results
                                           results.models <- list()
                                           
                                           #############################
                                           
                                           ### Tissue offset model ###
                                           glmmTMB.m1 <- glmmTMB(counts ~  nf + sex +  time * copd + (1|subject), 
                                                           
                                                                 offset = tissue, 
                                                                 data = dat,
                                                                 family = "nbinom1")
                                           
                                       
                                        
                                           pos.def.m1 <- TRUE
                                           
                                           # Check if models has positive definite Hessian matrix 
                                           if(!glmmTMB.m1$sdr$pdHess) {
                                                   
                                                   # If not, simplify fit to Poisson family...        
                                                   glmmTMB.m1 <- glmmTMB(counts ~  nf +  time * copd + (1|subject), 
                                                                         offset = tissue, 
                                                                         data = dat,
                                                                         family = "poisson")
    
                                                   pos.def.m1 <- FALSE
                                                   
                                           }
                                           
                           
                                        
                                        #  anova.table <-  data.frame(car::Anova(glmmTMB.m1))
                                        #   
                                     #
                                        #  global.p.m1 <- anova.table[4, 3]
                                           
                                           
                                           # Cui et al style model 
                                
                                          glmmTMB.m2 <- glmmTMB(counts ~  nf + sex + time * copd + (1|subject), 
                                                      
                                                                data = dat,
                                                                family = "nbinom1")
                                         
                                       
                                           pos.def.m2 <- TRUE
                                     
                                          # Check if models has positive definite Hessian matrix 
                                          if(!glmmTMB.m2$sdr$pdHess) {
                                                  
                                                  # If not, simplify fit to Poisson family...        
                                                  
                                                  
                                                  pos.def.m2 <- FALSE
                                                  
                                                  
                                                  glmmTMB.m2 <- glmmTMB(counts ~  nf +  time * copd + (1|subject), 
                                                                     
                                                                        data = dat,
                                                                        family = "poisson")
                                                  
                                                
                                          }
                                          
                          
                                          
                                         #  anova.table <-  data.frame(car::Anova(glmmTMB.m2))
                                         #  
                                         #  
                                         #  global.p.m2 <- anova.table[4, 3]
                                         #  
                                 
                                           
                                           
                                           if(diagnostics) {
                                                   
                                         
                                                   res.glmmTMB.m1   <-  simulateResiduals(glmmTMB.m1   , n = 1000)
                                           
                                                   res.glmmTMB.m2   <-  simulateResiduals(glmmTMB.m2   , n = 1000)
                                                   
                                             
                                                   disp.m1    <- testDispersion( res.glmmTMB.m1  , plot = FALSE)
                                                   disp.m2    <- testDispersion( res.glmmTMB.m2  , plot = FALSE)
                                                   
                                                  
                                                   uniform.m1 <-   testUniformity(res.glmmTMB.m1  , plot = FALSE)
                                                   uniform.m2 <-   testUniformity(res.glmmTMB.m2  , plot = FALSE)
                                                   
                                                   
                                                   
                                           } else {
                                                   
                                                   disp.m1   <- list(statistic = NA, p.value = NA)  
                                                   disp.m2   <- list(statistic = NA, p.value = NA) 
                                                   
                                                   
                                                   uniform.m1 <- list(statistic = NA, p.value = NA) 
                                                   uniform.m2 <- list(statistic = NA, p.value = NA) 
           
                                           }
                                           
                                           
                                           glmmTMB.results  <- rbind(data.frame(summary(glmmTMB.m1)$coef$cond) %>%
                                                   mutate(coef = rownames(.), 
                                                          model = "tissue_offset", 
                                                          interaction = TRUE, 
                                                          gene = genes[i], 
                                                          global.p = NA,
                                                          theta = sigma(glmmTMB.m1), 
                                                          dispersion.statistic = disp.m1$statistic, 
                                                          dispersion.p = disp.m1$p.value, 
                                                          uniform.statistic = uniform.m1$statistic,
                                                          uniform.p = uniform.m1$p.value, 
                                                          positive.def = pos.def.m1),
                                          
                                                   data.frame(summary(glmmTMB.m2)$coef$cond) %>%
                                                                              mutate(coef = rownames(.), 
                                                                                     model = "nf", 
                                                                                     interaction = TRUE, 
                                                                                     gene = genes[i], 
                                                                                     global.p = NA,
                                                                                     theta = sigma(glmmTMB.m2), 
                                                                                     dispersion.statistic = disp.m2$statistic, 
                                                                                     dispersion.p = disp.m2$p.value, 
                                                                                     uniform.statistic = uniform.m2$statistic,
                                                                                     uniform.p = uniform.m2$p.value,
                                                                                  
                                                                                     positive.def = pos.def.m2)) %>%
                                                   dplyr::select(gene, 
                                                                 model,
                                                                 interaction,
                                                                 coef, 
                                                                 estimate = Estimate, 
                                                                 se = Std..Error, 
                                                                 z.val = z.value, 
                                                                 p.val = Pr...z.., 
                                                                 global.p = global.p,
                                                                 theta,
                                                                 dispersion.statistic,
                                                                 dispersion.p,        
                                                                 uniform.statistic,   
                                                                 uniform.p, 
                                                                 positive.def) 
                                           
                                           
                                           results <- glmmTMB.results
                                           
                                           
                                           
                                           # return results
                                           results
                                           
                                   },
                                   error = function(e){
                                           message('** ERR at ', Sys.time(), " **")
                                           print(e)
                                           
                                           
                                   })
                           
                           
                           
                   }

close(pb)
stopCluster(cl)

results2 <- bind_rows(results)

### This has been saved!  
saveRDS(results2, file = "./data/derivedData/mixed_model/mixed_full_copd.RDS")
#rm(results)
#rm(results2)

