########### Cerno Gene ontology analysis ############

# The non-directinal rank based test of enrichment is performed with 
# the tmod package and gene ontologies. 
# A modified MSD metric is used as the gene level metric. 





source("./R/lib_fun.R")
source("./R/model_coefs.R")

unique(all_coefs$coef)


ac_d <-  all_coefs %>%
  filter(coef %in% c("copdcopd", "timeW3:copdcopd", "timePostExc:copdcopd")) %>% 
        mutate(ciu = estimate + qnorm(0.975) * se, 
               cil = estimate - qnorm(0.975) * se, 
               msd = if_else(estimate > 0, cil, -ciu), 
               log2fc = estimate/log(2), 
               model_coef = paste0(model, "_", coef)) %>%  # msd calculation Wald 95% CI
        
        # A modified msd: when msd < 0, genes are ranked based on unadjusted p-values
        mutate(msd.mod = if_else(msd > 0, msd, -p.val)) %>%
        
        print()



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




hallmark$MODULES$ID <- hallmark$MODULES$Title
kegg$MODULES$ID <- kegg$MODULES$Title

msig.bp$MODULES$ID <- msig.bp$MODULES$Title
msig.cc$MODULES$ID <- msig.cc$MODULES$Title
msig.mf$MODULES$ID <- msig.mf$MODULES$Title
 

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

temp.hallmark <- hallmark$MODULES %>%
  inner_join(count_gene_cat(hallmark$MODULES2GENES, symbols$SYMBOL))

temp.kegg <- kegg$MODULES %>%
  inner_join(count_gene_cat(kegg$MODULES2GENES, symbols$SYMBOL))

temp.bp <- msig.bp$MODULES %>%
        inner_join(count_gene_cat(msig.bp$MODULES2GENES, symbols$SYMBOL))

temp.cc <- msig.cc$MODULES %>%
        inner_join(count_gene_cat(msig.cc$MODULES2GENES, symbols$SYMBOL))

temp.mf <- msig.mf$MODULES %>%
        inner_join(count_gene_cat(msig.mf$MODULES2GENES, symbols$SYMBOL))


# Coefs for loop     
coefs <-  unique(ac_d$model_coef)     

# Store results in lists
cerno_results <- list()
u_results <- list()
z_results <- list()


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

            
        # Remove large and small gene set sizes based on set size after subseting
        
        sel.bp <-  temp.bp$set.size.present <= 500 & temp.bp$set.size.present >= 25
        sel.cc <-  temp.cc$set.size.present <= 500 & temp.cc$set.size.present >= 25
        sel.mf <-  temp.mf$set.size.present <= 500 & temp.mf$set.size.present >= 25
        
        
        # Perform cerno algorithm 
        res.hallmark <- tmodCERNOtest(gl_symb$SYMBOL, mset = hallmark, qval = 1)
        res.kegg <- tmodCERNOtest(gl_symb$SYMBOL, mset = kegg, qval = 1)
        res.bp <- tmodCERNOtest(gl_symb$SYMBOL, mset = msig.bp[sel.bp], qval = 1)
        res.cc <- tmodCERNOtest(gl_symb$SYMBOL, mset = msig.cc[sel.cc], qval = 1)
        res.mf <- tmodCERNOtest(gl_symb$SYMBOL, mset = msig.mf[sel.mf], qval = 1)
        
      
        
        cerno <- res.bp %>%
                tibble::rownames_to_column(var = "id") %>%
                dplyr::select(-id) %>%
                mutate(go.cat = "bp") %>%
                rbind(res.cc %>%
                              tibble::rownames_to_column(var = "id") %>%
                              dplyr::select(-id) %>%
                              mutate(go.cat = "cc")) %>%
                rbind(res.mf %>%
                              tibble::rownames_to_column(var = "id") %>%
                              dplyr::select(-id) %>%
                              mutate(go.cat = "mf")) %>%
          rbind(res.hallmark %>%
                  tibble::rownames_to_column(var = "id") %>%
                  dplyr::select(-id) %>%
                  mutate(go.cat = "hallmark")) %>%
          rbind(res.kegg %>%
                  tibble::rownames_to_column(var = "id") %>%
                  dplyr::select(-id) %>%
                  mutate(go.cat = "kegg")) %>%
               
           mutate(coef = coefs[i])
        
      
        
        cerno_results[[i]] <- cerno
        
        message("Iteration ", i , " of ", length(coefs))
        
}

### Save results ###


cerno_df <- bind_rows(cerno_results)

saveRDS(cerno_df, "./data/derivedData/cerno.RDS")


