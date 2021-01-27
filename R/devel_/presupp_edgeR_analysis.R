##-------------------------------------
## presupp_edgeR_analysis.R
##
## Title:
## Purpose: 
## Author:
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



### Compare baseline estrimates between edgeR and mixed model (presupp vs. preexc)

# Load edgeR
results_edgeR <- readRDS(file = "./data/derivedData/edgeR/presupp_results.RDS")


# named vectors of gene-level statistics
gene_name_convert <- clusterProfiler::bitr(unique(results_edgeR$gene), fromType = "ENSEMBL",
                                           toType = c( "ENTREZID", "SYMBOL"),
                                           OrgDb = org.Hs.eg.db)


er_results <- results_edgeR %>%
  inner_join(gene_name_convert %>%
               dplyr::select(gene = ENSEMBL, SYMBOL, ENTREZID)) %>%
  print()





### GSEA with Hallmark 



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




ac_d <-  er_results %>%

  mutate(msd = if_else(logFC > 0, cil, -ciu), 
         log2fc = logFC) %>%  # msd calculation Wald 95% CI
  
  # A modified msd: when msd < 0, genes are ranked based on unadjusted p-values
  mutate(msd.mod = if_else(msd > 0, msd, -PValue)) %>%
  # Pi-value, combining log fold change and p-values
  # from: Xiao et al A novel significance score for gene selection and ranking
  mutate(pi = log2fc * -log10(FDR)) %>%
  
  print()







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


# Unique symbol identifiers from the data set
symbols <- data.frame(SYMBOL = unique(ac_d$SYMBOL))


# Create temporary data frame with set sizes and set size that are present
hallmark_gs <- hallmark$MODULES %>%
  inner_join(count_gene_cat(hallmark$MODULES2GENES, symbols$SYMBOL))


gl <- ac_d %>%
  ungroup() %>%
  arrange(-msd.mod) %>%
  pull(gene)


# Convert gene id to symbols for use in 
gl_symb <- bitr(gl, fromType = "ENSEMBL",
                toType = c( "ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)



# Perform cerno algorithm 
res.hallmark.cc <- tmodCERNOtest(gl_symb$SYMBOL, mset = hallmark, qval = 1)



cerno <- res.hallmark.cc %>%
  tibble::rownames_to_column(var = "id") %>%
  dplyr::select(-id) %>%
  mutate(go.cat = "hallmark") 



### GSEA

translated_data <- ac_d %>%
  ungroup() %>%
  inner_join(gene_name_convert %>% 
               mutate(gene = ENSEMBL)) %>%
  dplyr::select(gene,
                ENTREZID, 
                SYMBOL, 
      
                estimate = logFC,
                log2fc, 
                pi) %>%
  group_by(SYMBOL) %>%
  summarise(estimate = mean(pi)) 


est <- translated_data %>%
  pull(estimate)

names(est) <- translated_data %>%
  pull(SYMBOL)



gene_set_list_hallmark <- hallmark$MODULES2GENES
## Gene list 1: Hallmark
fgseaRes.hallmark <- fgsea(pathways = gene_set_list_hallmark, 
                           stats = est, 
                           minSize = 0, 
                           maxSize = 500,
                           eps = 0,
                           nproc = 8) 


fgsea_full <- fgseaRes.hallmark %>%
  data.frame() %>%
  mutate(go.cat = "hallmark") 
  


gsea <- fgsea_full %>%
  dplyr::select(ID = pathway,
                gsea.padj = padj,
                gsea.p = pval, 
                gsea.es = ES, 
                gsea.nes = NES) 



# full cerno
cerno <- cerno %>%
  dplyr::select(ID = Title,
                cerno, 
                N1, 
                cerno.p = P.Value,
                cerno.padj = adj.P.Val, 
                go.cat) 



## Set significance level for gene set analysis
sig.lev <- 0.05

# Top fgsea gene categories
comb <- gsea %>%
  inner_join(cerno) %>%
  mutate(sig = if_else(gsea.padj < sig.lev & cerno.padj < sig.lev, "con", 
                       if_else(gsea.padj < sig.lev &  cerno.padj > sig.lev, "fgsea", 
                               if_else(gsea.padj > sig.lev &  cerno.padj < sig.lev, "cerno", "non")))) %>% 
  filter(sig != "non") %>%
  print()
  
# Change coef into two variables 




