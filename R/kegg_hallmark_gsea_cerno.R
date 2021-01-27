##-------------------------------------
##  kegg_hallmark_gsea_cerno.R
##
## Title: KEGG and hallmark gene set collections analyzed with GSEA and cerno
## Purpose: Perform gene set analysis using kegg and hallmark gene set collections
## Author: Daniel Hammarström
##
##
##
##-------------------------------------
## Notes: 
#
#
#
#
#
#
#
#
## ------------------------------------



source("./R/lib_fun.R")

# Download gene sets from msigdb --> 
# Hallmark and KEGG gene sets used for higher level categories to complement gene ontology. 



hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
        dplyr::select(ID  = gene_symbol, 
                      Title = gs_name) %>%
        data.frame() %>%
        makeTmodFromDataFrame(feature_col  = 1, 
                              module_col = 2) %>%
        print()


biocarta <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
        dplyr::select(ID  = gene_symbol, 
                      Title = gs_name) %>%
        data.frame() %>%
        makeTmodFromDataFrame(feature_col  = 1, 
                              module_col = 2) %>%
        print()




# Load all models and coefs 
train_mixed <- readRDS(file = "./data/derivedData/DE/mixed_separate.RDS")

presupp_mixed <- readRDS(file = "./data/derivedData/DE/mixed_presupp.RDS")


mixed_results_min <- rbind(train_mixed, presupp_mixed) %>%
        print()

# Load raw counts files. dge_sum is prepared in seq_prepare.R
dge <- readRDS(file = "./data/derivedData/dge_lists/dge_sum.RDS")




#### Gene ontology analysis ###################

## Create data set for rank based tests      

## Checking individual coefficients ##

all_coefs <- mixed_results_min %>%
        group_by(coef, interaction, model) %>%
        mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
        ungroup() %>%
  filter(uniform.p > 0.05) %>%
        dplyr::select(gene, model, interaction, coef, estimate, se, z.val, p.val, p.adj) %>%
        print()


ac_d <-  all_coefs %>%
        filter((coef %in% c("timeW3:treatvitd",
                                   "timePostExc:treatvitd",
                                   "timePreExc:treatvitd"))) %>% 
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





# Retreive saved data

msig.list <- readRDS(file = "./data/derivedData/msig_tmod_lists.RDS")

msig.bp <- msig.list$msig.bp
msig.cc <- msig.list$msig.cc
msig.mf <- msig.list$msig.mf

msig.bp$MODULES$ID <- msig.bp$MODULES$Title
msig.cc$MODULES$ID <- msig.cc$MODULES$Title
msig.mf$MODULES$ID <- msig.mf$MODULES$Title

rownames(msig.bp$MODULES) <- names(msig.bp$MODULES2GENES) 
rownames(msig.cc$MODULES) <- names(msig.cc$MODULES2GENES) 
rownames(msig.mf$MODULES) <- names(msig.mf$MODULES2GENES) 


# fgsea filters gene sets based on number of genes present in each category as 
# filtering criteria. The modules has a B variable containing the number of 
# genes in each category, some of these are missing from our data set. 
# A new B variable is created to fix this. 

# A list of symbols corresponding to our data set #
# Convert gene id to symbols for use in 

symbols <- bitr(unique(ac_d$gene), fromType = "ENSEMBL",
                toType = c( "ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

# The function counts the number of genes in each gene set after removing genes in gene sets 
# that are not present/detected in our data set.

count_gene_cat <- function(x, symbols){
        
        results <- data.frame(ID = names(x), 
                              set.size = rep(NA, length(x)), 
                              set.size.present = rep(NA, length(x)))
        
        for(i in 1:nrow(results)) {
                
                results[i, 2] <- length(x[[results[i,1]]])  
                results[i, 3] <- sum(x[[results[i,1]]] %in% symbols)    
                
        }
        return(results)
}  

# Create temporary data frame with set sizes and set size that are present
hallmark_gs <- hallmark$MODULES %>%
        inner_join(count_gene_cat(hallmark$MODULES2GENES, symbols$SYMBOL))

biocarta_gs <- biocarta$MODULES %>%
        inner_join(count_gene_cat(biocarta$MODULES2GENES, symbols$SYMBOL))


# Coefs for loop     
coefs <-  unique(ac_d$model_coef)     

# Store results in lists
cerno_results <- list()



for(i in 1:length(coefs)) {
        
        
        # Cerno algorithm 
        
        gl <- ac_d %>%
                ungroup() %>%
                filter(model_coef == coefs[i]) %>%
                arrange(-msd.mod) %>%
                pull(gene)
        
        
        # Convert gene id to symbols for use in 
        gl_symb <- bitr(gl, fromType = "ENSEMBL",
                        toType = c( "ENTREZID", "SYMBOL"),
                        OrgDb = org.Hs.eg.db)
        
        
        
        # Perform cerno algorithm 
        res.hallmark.cc <- tmodCERNOtest(gl_symb$SYMBOL, mset = hallmark, qval = 1)
        res.biocarta.cc <- tmodCERNOtest(gl_symb$SYMBOL, mset = biocarta, qval = 1)

        
        
        cerno <- res.hallmark.cc %>%
                tibble::rownames_to_column(var = "id") %>%
                dplyr::select(-id) %>%
                mutate(go.cat = "hallmark") %>%
                rbind(res.biocarta.cc %>%
                              tibble::rownames_to_column(var = "id") %>%
                              dplyr::select(-id) %>%
                              mutate(go.cat = "biocarta")) %>%
                mutate(coef = coefs[i])
        
    
        
        cerno_results[[i]] <- cerno
        
        
        message("Iteration ", i , " of ", length(coefs))
        
}

### Save results ###


cerno_df <- bind_rows(cerno_results)

saveRDS(list(cerno.test = cerno_df), 
        "./data/derivedData/cerno_hallmark_biocarta.RDS")


########### Gene set enrichment analysis ###########################


# Prepare gene sets ################
# Gene sets downloaded from MSigDG. 

#### Using msig list from cerno_analysis.R
msig.list <- readRDS(file = "./data/derivedData/msig_tmod_lists.RDS")

msig.bp <- msig.list$msig.bp
msig.cc <- msig.list$msig.cc
msig.mf <- msig.list$msig.mf

# Convert to fit the fgsea function/package
names(msig.bp$MODULES2GENES) <- msig.bp$MODULES$Title
names(msig.cc$MODULES2GENES) <- msig.cc$MODULES$Title
names(msig.mf$MODULES2GENES) <- msig.mf$MODULES$Title


gene_set_list_hallmark <- hallmark$MODULES2GENES
gene_set_list_biocarta <- biocarta$MODULES2GENES


## Combine all in a list
gene_lists <- list(hallmark = gene_set_list_hallmark, 
                   biocarta = gene_set_list_biocarta) 


# Convert 

# named vectors of gene-level statistics

gene_name_convert <- bitr(unique(ac_d$gene), fromType = "ENSEMBL",
                          toType = c( "ENTREZID", "SYMBOL"),
                          OrgDb = org.Hs.eg.db)


results_models_collapsed <- list()
results_models_full <- list()

coefs <- unique(ac_d$model_coef)


full_gsea <- list()
collapse_gsea <- list()

ac_d %>%
        filter(model == "nf", 
               coef %in% c("timeW3:treatvitd")) %>%
        inner_join(gene_name_convert %>% 
                           mutate(gene = ENSEMBL)) %>%
        filter(SYMBOL %in% gene_set_list_biocarta$BIOCARTA_VDR_PATHWAY) %>%
        print()



### For loop over gene sets using fgsea

for(i in 1:length(coefs)) { 
        
        #### Beginning of loop tasks                  
        
        
        translated_data <- ac_d %>%
                filter(model_coef == coefs[i]) %>%
                ungroup() %>%
                inner_join(gene_name_convert %>% 
                                   mutate(gene = ENSEMBL)) %>%
                dplyr::select(gene,
                              ENTREZID, 
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
        fgseaRes.hallmark <- fgsea(pathways =  gene_lists[[1]], 
                             stats = est, 
                             minSize = 0, 
                             maxSize = 500,
                             eps = 0,
                             nproc = 8) 
        

        
        ## Gene list 2: Biocarta
        fgseaRes.biocarta <- fgsea(pathways =  gene_lists[[2]], 
                             stats = est, 
                             minSize = 0, 
                             maxSize = 500,
                             eps = 0,
                             nproc = 8) 
        
  
        
        
        
  
        ### Combine results 
        
        fgsea_full <- fgseaRes.hallmark %>%
                data.frame() %>%
                mutate(go.cat = "hallmark") %>%
                rbind(fgseaRes.biocarta %>%
                              data.frame() %>%
                              mutate(go.cat = "biocarta")) %>%
                mutate(coef = coefs[i]) 
        
        

        
        
        full_gsea[[i]] <- fgsea_full
        

        message(paste0("Model coefficient: ", i, "/", length(coefs) ))
        
}

# Combine results 

full_gsea_df <- bind_rows(full_gsea)


saveRDS(full_gsea_df, file = "./data/derivedData/fgsea_full_hallmark_biocarta.RDS")


full_gsea_df %>%
        dplyr::select(pathway:size, go.cat, coef) %>%
        filter(go.cat == "biocarta") %>%
        arrange(padj) %>%
        print()
        



