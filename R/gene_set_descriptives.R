

### Descriptive data for gene sets ####

# This script produces a data frame with gene ontology stats are expressed over
# model/coefs

# E.g. n genes with msd > 0, number of genes in set, n genes with log2fc > | < 0. 



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

gene_name_convert <- clusterProfiler::bitr(unique(ac_d$gene), fromType = "ENSEMBL",
                          toType = c( "ENTREZID", "SYMBOL"),
                          OrgDb = org.Hs.eg.db)

# df with symbols
symbols_results <- ac_d %>%
        inner_join(gene_name_convert %>%
                           mutate(gene  = ENSEMBL) %>%
                           dplyr::select(gene, SYMBOL)) %>%
        mutate(coefs = paste0(model, "_", coef)) %>% 
        print()






##### Load fgsea results to retrieve leading edge genes 

# full gsea
full_gsea_df <- readRDS(file = "./data/derivedData/fgsea_full.RDS")

gsea <- full_gsea_df %>%
        dplyr::select(ID = pathway,
                      leadingEdge,
                      coef) 


# Filter gene sets 

gene_lists$bp <-       gene_lists$bp[names(gene_lists$bp) %in% gsea$ID]
gene_lists$mf <-       gene_lists$mf[names(gene_lists$mf) %in% gsea$ID]
gene_lists$cc <-       gene_lists$cc[names(gene_lists$cc) %in% gsea$ID]
gene_lists$hallmark <- gene_lists$hallmark[names(gene_lists$hallmark) %in% gsea$ID]
gene_lists$kegg     <- gene_lists$kegg[names(gene_lists$kegg) %in% gsea$ID]


names(gene_lists)
names(gl) == 'Go inner mitochondrial membrane protein complex'
i <- 2
j <- 153
k <- 9
# for every gene set list (gene_lists i), calculate summary statistics in
# every gene set (j) for every coefficient/model (k). Compile in data frames.

results <- list()

for(i in 1:length(gene_lists))  {
        
     
        gl <- gene_lists[[i]]

        gene_list_results <- list()
        
        for(j in 1:length(gl)){
               
                gl.name <- names(gl)[j]
                
                coefs <- unique(symbols_results$coefs)
             
                coefs_results <- list() 
                
         
                
                for(k in 1:length(coefs)){
                        
                        # Get leading edge genes from GSEA
                        # These genes contribute to the GSEA score
                        # See ?fgsea
                        leading_edge <- gsea %>%
                                filter(ID == gl.name, 
                                       coef == coefs[k]) %>%
                                pull(leadingEdge) %>%
                                unlist()
                        
                        # Explicit selection of coefficient!
                        select.coef <- coefs[k]
                        
                        coefs_results[[k]] <- symbols_results %>%
                                        
                               filter(coefs == select.coef) %>%
                              
                                filter(SYMBOL %in% gl[[gl.name]]) %>%
                             
                                mutate(msd.pos = if_else(msd > 0, 1, 0), 
                                       in_LE = if_else(SYMBOL %in% leading_edge, 1, 0), 
                                       msd.pos_inLE = if_else(msd.pos == 1 & in_LE == 1, 1, 0), 
                                       fc.in_LE = if_else(SYMBOL %in% leading_edge, log2fc, NULL)) %>%
                                        
                                
                                dplyr::select(gene, log2fc, msd, SYMBOL, msd.pos, in_LE, msd.pos_inLE, fc.in_LE) %>%
                 
                                        
                                group_by() %>%
                                summarise(msd.pos = sum(msd.pos, na.rm = TRUE), # n genes of gene set with positive msd
                                          in_LE = sum(in_LE),  # n genes in gene set in leading edge
                                          msd.pos_inLE = sum(msd.pos_inLE, na.rm = TRUE), # n genes in leading edge with msd > 0 
                                          n.total = n(), # n total genes in identified from gene set 
                                          
                                          fc.LE.min = min(fc.in_LE, na.rm = TRUE), 
                                          fc.LE.max = max(fc.in_LE, na.rm = TRUE), 
                                          fc.LE = mean(fc.in_LE, na.rm = TRUE))  %>% 
            
                                mutate(set.size = length(gl[[j]]), # n genes in raw gene set
                                       ID = gl.name, 
                                       coefs = select.coef, 
                                       go.cat = names(gene_lists)[i]) 
                            
                }
   
      gene_list_results[[j]] <-  bind_rows(coefs_results)    
      print(paste0("Gene set ",j , " of ", length(gl), " in gene list ", i))  
                 
        }
        results[[i]] <- bind_rows(gene_list_results)
}


results2 <- bind_rows(results)


saveRDS(results2, file = "./data/derivedData/geneset_descriptives.RDS")




















