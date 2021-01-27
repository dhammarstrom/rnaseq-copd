##-------------------------------------
## mitocarta_analysis.R
##
## Title: 
## Purpose: 
## Author:
##
##
##
##-------------------------------------
## Notes:
# A sub analysis is perfomed on mitochondrial genes for the puspose of
# exploring differences between COPD/Healthy. 
#
#
#
#
#
#
## ------------------------------------

library(mitocarta)
source("./R/lib_fun.R")




# Genes with known mitochondrial localization in skeletal muscle are
# retrieved from mitocharta (v. 2.0, possibly update to 3.0 later)


mito_skm <- read_excel("./data/annotations/Human.MitoCarta3.0.xls", 
                        sheet = "B Human All Genes") %>%
  dplyr::select(HumanGeneID, Symbol, 
                skeletalmuscle_total_peak_intensity_log10) %>%
  filter(skeletalmuscle_total_peak_intensity_log10 > 0) %>%

  distinct(Symbol) %>%
  pull(Symbol)


# mitocarta 3.0

mito_pathways <- read_excel("./data/annotations/Human.MitoCarta3.0.xls", 
                         sheet = "C MitoPathways") %>%
  dplyr::select("MitoPathway", 
                "Genes") %>%
  data.frame()


results <- list()

for(i in 1:nrow(mito_pathways)) {
  
  results[[i]] <- data.frame(ID = unlist(strsplit(mito_pathways[i,2], split = ", ")), 
             Title = mito_pathways[i,1])

}

mito_pathways_df <- bind_rows(results)



mito_module <- mito_pathways_df %>%
  filter(ID %in% mito_skm) %>%
  makeTmodFromDataFrame(feature_col  = 1, 
                        module_col = 2) %>%
  print()



# Get gene-wise estimates from mixed models with COPD condition interaction with training

## Load DE model data

mixed_dat <- readRDS(file = "./data/derivedData/mixed_model/mixed_full_copd.RDS")

## Convert gene ids to symbol

gene_name_convert <- clusterProfiler::bitr(unique(mixed_dat$gene), fromType = "ENSEMBL",
                                           toType = c( "ENTREZID", "SYMBOL"),
                                           OrgDb = org.Hs.eg.db) %>%
  dplyr::select(gene = ENSEMBL, SYMBOL)


mixed_dat %>%
  filter(gene == "ENSG00000184470")

all_coefs_df <- mixed_dat %>%
  
  filter(uniform.p > 0.05) %>%
  
  inner_join(gene_name_convert) %>%
  
  filter(SYMBOL %in% mito_skm) %>%

  group_by(coef, model, interaction) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr"), 
         log2fc = estimate/log(2)) %>%
  ungroup() %>%
 filter(model == "nf", 
         coef %in% c("copdcopd", 
                     "timePostExc:copdcopd", 
                     "timeW3:copdcopd")) %>%

  mutate(ciu = estimate + qnorm(0.975) * se, 
         cil = estimate - qnorm(0.975) * se, 
         msd = if_else(estimate > 0, cil, -ciu), 
         log2fc = estimate/log(2), 
         model_coef = paste0(model, "_", coef)) %>%  # msd calculation Wald 95% CI
  
  # A modified msd: when msd < 0, genes are ranked based on unadjusted p-values
  mutate(msd.mod = if_else(msd > 0, msd, -p.val)) %>%
  # Pi-value, combining log fold change and p-values
  # from: Xiao et al A novel significance score for gene selection and ranking
  mutate(pi = log2fc * -log10(p.val)) %>%
  
  print()
  
  
coefs <- unique(all_coefs_df$coef)

cerno_results <- list()
full_gsea <- list()
  
for(i in 1:length(coefs)) {
    
    
    # Cerno algorithm 
    
    gl <- all_coefs_df %>%
      ungroup() %>%
      filter(coef == coefs[i]) %>%
      arrange(-msd.mod) %>%
      pull(gene)
    
    
    # Convert gene id to symbols for use in 
    gl_symb <- bitr(gl, fromType = "ENSEMBL",
                    toType = c( "ENTREZID", "SYMBOL"),
                    OrgDb = org.Hs.eg.db)
    
    
    
    # Perform cerno algorithm 
    res.mito.cc <- tmodCERNOtest(gl_symb$SYMBOL, mset = mito_module, qval = 1)

    
    
    
    cerno <- res.mito.cc %>%
      tibble::rownames_to_column(var = "id") %>%
      dplyr::select(-id) %>%
      mutate(go.cat = "mito") %>%

      mutate(coef = coefs[i])
    
    
    
    cerno_results[[i]] <- cerno
    
    
    #### FGSEA 
    
    translated_data <- all_coefs_df %>%
      filter(coef == coefs[i]) %>%
      ungroup() %>%
      inner_join(gene_name_convert) %>%
      dplyr::select(gene,
                  
                    SYMBOL, 
                    z.val, 
                    estimate,
                    log2fc, 
                    pi) %>%
      group_by(SYMBOL) %>%
      summarise(estimate = mean(pi)) 
    
    
    est <- translated_data %>%
      pull(estimate)
    
    names(est) <- translated_data %>%
      pull(SYMBOL)
    
    ## Gene list 1: Hallmark
    fgseaRes.mito <- fgsea(pathways =  mito_module$MODULES2GENES, 
                               stats = est, 
                               minSize = 0, 
                               maxSize = 1000,
                               eps = 0,
                               nproc = 8) 
    
 
    ### Combine results 
    
    fgsea_full <- fgseaRes.mito %>%
      data.frame() %>%
      mutate(go.cat = "mito") %>%

      mutate(coef = coefs[i]) 
    
    
    
    
    
    full_gsea[[i]] <- fgsea_full
    
    
    message(paste0("Model coefficient: ", i, "/", length(coefs) ))
    
    
  }

### Save results ###


cerno_df <- bind_rows(cerno_results)
fgsea_df <- bind_rows(full_gsea)


saveRDS(list(cerno.mito = cerno_df, 
             fgsea.mito = fgsea_df), 
        "./data/derivedData/cerno_fgsea_mitocarta.RDS")

  
  
  
################################# Mitocarta per condition models ################

# Using only mitocarta genes, posthoc tests are computed per group from baseline


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

mito_skm

gene_name_convert <- clusterProfiler::bitr(mito_skm, fromType = "SYMBOL",
                                           toType = c("ENSEMBL",
                                                      "ENTREZID"),
                                           OrgDb = org.Hs.eg.db) %>%
  dplyr::select(gene = ENSEMBL, SYMBOL)



# Set what dge list to use
selected.dge <- dge

# A list of genes (available in the dge-list)
genes <- rownames(selected.dge)[rownames(selected.dge)  %in% gene_name_convert$gene]

# set up parallel processing using the foreach package 
cores <- as.integer(round(detectCores() * 1, 0))
cl <- makeCluster(cores[1]) # not to overload cpu #### Change this if on server!
registerDoSNOW(cl)

# Progress bar (should work on linux also)
iterations <- length(genes) 
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)



# foreach loop 
results <- foreach(i = 1:length(genes), 
                   .packages = c("glmmTMB", "dplyr", "DHARMa", "car", "emmeans"),  # Packages used in fitting/data wrangling
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
                         
                        
                         # Create a control group
                         
                 
                         ## Extract gene counts and put together in the data frame
               
                         
                         
                         
                          ## Extract gene counts and put together in the data frame
                         dat <- data.frame(counts = selected.dge$counts[genes[i],], 
                                           sample = colnames(selected.dge$counts)) %>%
                           inner_join(dat) %>%
                           inner_join(dat_yield) %>%
                           inner_join(subject_sex) %>%
                           filter(!(time %in% c("PreSupp"))) %>%
                           
                           mutate(counts = as.integer(round(counts, 0)), 
                                  subject = factor(subject), 
                                  treat = factor(treat, levels = c("placebo", "vitd")), 
                                  time = factor(time, levels = c("PreExc",
                                                                 "W3",
                                                                 "PostExc")),
                                  leg = paste0(leg, subject),
                                  eff.lib = lib.size * norm.factors, 
                                  tissue = log(1000/yield))
                         
                         
                        
                         
                         # A list of results
                         results.models <- list()
                         
                         #############################
                         
                         ### Tissue offset model ###
                         glmmTMB.m1 <- glmmTMB(counts ~  nf + sex +  time + (1|subject), 
                                               
                                               offset = tissue, 
                                               data = dat,
                                               family = "nbinom1")
                         
                         
                       
                         pos.def.m1 <- TRUE
                         
                         # Check if models has positive definite Hessian matrix 
                         if(!glmmTMB.m1$sdr$pdHess) {
                           
                           # If not, simplify fit to Poisson family...        
                           glmmTMB.m1 <- glmmTMB(counts ~  nf +  time + (1|subject), 
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
                         
                         glmmTMB.m2 <- glmmTMB(counts ~  nf + sex + time  + (1|subject), 
                                               
                                               data = dat,
                                               family = "nbinom1")
                         
                         
                         pos.def.m2 <- TRUE
                         
                         # Check if models has positive definite Hessian matrix 
                         if(!glmmTMB.m2$sdr$pdHess) {
                           
                           # If not, simplify fit to Poisson family...        
                           
                           
                           pos.def.m2 <- FALSE
                           
                           
                           glmmTMB.m2 <- glmmTMB(counts ~  nf +  time  + (1|subject), 
                                                 
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


  
  
  