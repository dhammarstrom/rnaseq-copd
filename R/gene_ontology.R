########### GENE ONTOLOGY ANALYSIS #############

# Using fgsea making gene set enrichment analysis based on 
# pi is the gene level metric (Xiao et al A novel significance score for gene selection and ranking).
# Interaction terms are used to make comparisons.


# Prepare data set 
source("./R/lib_fun.R")
source("./R/model_coefs.R")


ac_d <-  all_coefs %>%
    filter(coef %in% c("copdcopd", "timeW3:copdcopd", "timePostExc:copdcopd")) %>% 
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




# Prepare gene sets ################
# Gene sets downloaded from MSigDG. 
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(ID  = gene_symbol, 
                  Title = gs_name) %>%
    data.frame() %>%
    makeTmodFromDataFrame(feature_col  = 1, 
                          module_col = 2) %>%
    print()


kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
    dplyr::select(ID  = gene_symbol, 
                  Title = gs_name) %>%
    data.frame() %>%
    makeTmodFromDataFrame(feature_col  = 1, 
                          module_col = 2) %>%
    print()

msig.bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
    dplyr::select(ID  = gene_symbol, 
                  Title = gs_name) %>%
    data.frame() %>%
    makeTmodFromDataFrame(feature_col  = 1, 
                          module_col = 2) %>%
    print()

msig.cc <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC") %>%
    dplyr::select(ID  = gene_symbol, 
                  Title = gs_name) %>%
    data.frame() %>%
    makeTmodFromDataFrame(feature_col  = 1, 
                          module_col = 2) %>%
    print()

msig.mf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
    dplyr::select(ID  = gene_symbol, 
                  Title = gs_name) %>%
    data.frame() %>%
    makeTmodFromDataFrame(feature_col  = 1, 
                          module_col = 2) %>%
    print()




hallmark$MODULES$ID<- hallmark$MODULES$Title
kegg$MODULES$ID    <- kegg$MODULES$Title
msig.bp$MODULES$ID <- msig.bp$MODULES$Title
msig.cc$MODULES$ID <- msig.cc$MODULES$Title
msig.mf$MODULES$ID <- msig.mf$MODULES$Title

## Combine all in a list
gene_lists <- list(bp =msig.bp$MODULES2GENES, 
                   cc =msig.cc$MODULES2GENES, 
                   mf =msig.mf$MODULES2GENES, 
                   hallmark = hallmark$MODULES2GENES, 
                   kegg = kegg$MODULES2GENES)


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
        
        ## Gene list 1: BP
        fgseaRes.bp <- fgsea(pathways =  gene_lists[[1]], 
                             stats = est, 
                             minSize = 25, 
                             maxSize = 500,
                             eps = 0,
                             nproc = 8) 
    

        
        ## Gene list 2: CC
        fgseaRes.cc <- fgsea(pathways =  gene_lists[[2]], 
                             stats = est, 
                             minSize = 25, 
                             maxSize = 500,
                             eps = 0,
                             nproc = 8) 
        
        ## Gene list 3: mf
        fgseaRes.mf <- fgsea(pathways =  gene_lists[[3]], 
                             stats = est, 
                             minSize = 25, 
                             maxSize = 500,
                             eps = 0,
                             nproc = 8) 
        
        ## Gene list 4: hallmark
        fgseaRes.hallmark <- fgsea(pathways =  gene_lists[[4]], 
                             stats = est, 
                             minSize = 0, 
                             maxSize = 500,
                             eps = 0,
                             nproc = 8) 
        
        ## Gene list 5: kegg
        fgseaRes.kegg <- fgsea(pathways =  gene_lists[[5]], 
                             stats = est, 
                             minSize = 0, 
                             maxSize = 500,
                             eps = 0,
                             nproc = 8) 
        
        
        
        ## Collapse fgsea results and store results from fgsea only
        
        collapsedPathways.bp <- collapsePathways(fgseaRes.bp[fgseaRes.bp$padj < 0.01, ],
                                                 pathways = gene_lists[[1]], 
                                                 stats = est)
        
        fgseaNonRedundant.bp <- fgseaRes.bp[fgseaRes.bp$pathway %in% collapsedPathways.bp$mainPathways, ]
        
        
        collapsedPathways.cc <- collapsePathways(fgseaRes.cc[fgseaRes.cc$padj < 0.01, ],
                                                 pathways = gene_lists[[2]], 
                                                 stats = est)
        
        fgseaNonRedundant.cc <- fgseaRes.cc[fgseaRes.cc$pathway %in% collapsedPathways.cc$mainPathways, ]
        
        
        collapsedPathways.mf <- collapsePathways(fgseaRes.mf[fgseaRes.mf$padj < 0.01, ],
                                                 pathways = gene_lists[[3]], 
                                                 stats = est)
        
        fgseaNonRedundant.mf <- fgseaRes.mf[fgseaRes.mf$pathway %in% collapsedPathways.mf$mainPathways, ]
        
        
        
        collapsedPathways.hallmark <- collapsePathways(fgseaRes.hallmark[fgseaRes.hallmark$padj < 0.01, ],
                                                 pathways = gene_lists[[4]], 
                                                 stats = est)
        
        fgseaNonRedundant.hallmark <- fgseaRes.hallmark[fgseaRes.hallmark$pathway %in% collapsedPathways.hallmark$mainPathways, ]
        
        
        collapsedPathways.kegg <- collapsePathways(fgseaRes.kegg[fgseaRes.kegg$padj < 0.01, ],
                                                       pathways = gene_lists[[5]], 
                                                       stats = est)
        
        fgseaNonRedundant.kegg <- fgseaRes.kegg[fgseaRes.kegg$pathway %in% collapsedPathways.kegg$mainPathways, ]
        
        
        
        
        
        ### Combine results 
        
        fgsea_full <- fgseaRes.bp %>%
                data.frame() %>%
                mutate(go.cat = "bp") %>%
                rbind(fgseaRes.cc %>%
                              data.frame() %>%
                              mutate(go.cat = "cc")) %>%
                rbind(fgseaRes.mf %>%
                              data.frame() %>%
                              mutate(go.cat = "mf")) %>%
            
            rbind(fgseaRes.hallmark %>%
                      data.frame() %>%
                      mutate(go.cat = "hallmark")) %>%
            
            rbind(fgseaRes.kegg %>%
                      data.frame() %>%
                      mutate(go.cat = "kegg")) %>%
            
            
            
                mutate(coef = coefs[i]) 
        
        

        fgsea_collapsed <- fgseaNonRedundant.bp %>%
                data.frame() %>%
                mutate(go.cat = "bp") %>%
                rbind(fgseaNonRedundant.cc %>%
                              data.frame() %>%
                              mutate(go.cat = "cc")) %>%
                rbind(fgseaNonRedundant.mf %>%
                              data.frame() %>%
                              mutate(go.cat = "mf")) %>%
            
            rbind(fgseaNonRedundant.hallmark %>%
                      data.frame() %>%
                      mutate(go.cat = "hallmark")) %>%
            
            rbind(fgseaNonRedundant.kegg %>%
                      data.frame() %>%
                      mutate(go.cat = "kegg")) %>%

                mutate(coef = coefs[i]) 
        
        
        
        
        full_gsea[[i]] <- fgsea_full
        
        collapse_gsea[[i]] <- fgsea_collapsed
        
        
        message(paste0("Model coefficient: ", i, "/", length(coefs) ))
        
}

# Combine results 

full_gsea_df <- bind_rows(full_gsea)
collapse_gsea <- bind_rows(collapse_gsea)


saveRDS(collapse_gsea, file = "./data/derivedData/fgsea_collapse.RDS")
saveRDS(full_gsea_df, file = "./data/derivedData/fgsea_full.RDS")



