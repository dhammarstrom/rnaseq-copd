##-------------------------------------
## figure-2.R
##
## Title: Mitocarta figure
## Purpose: Display within-group changes in mitocarat genes 
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

source("./R/lib_fun.R")
library(UpSetR)





gsea <- readRDS("./data/derivedData/cerno_fgsea_mitocarta.RDS")


gsea$cerno.mito %>%
  data.frame()%>%
  filter(coef == "timePostExc:copdcopd") %>%
  arrange(adj.P.Val) %>%

  print()


gsea$fgsea.mito %>%
  filter(coef == "timeW3:copdcopd") %>% 
  dplyr::select(-leadingEdge) %>%
  arrange(padj) %>%
  print()





### This has been saved!  



mm_results <- readRDS(
  file = "./data/derivedData/mixed_model/mixed_full_copd_mitocarta.RDS")


sig_genes <- mm_results$post_hoc %>%
  filter(contrast != "PostExc - W3") %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  filter(p.adj < 0.05, 
         estimate > 0) %>%
  mutate(coef = paste0(gsub(" - PreExc", "", contrast), 
                      "_", copd)) %>%
  
  
  print()


coefs <- unique(sig_genes$coef)
results <- list()


for(i in 1:length(coefs)) {
  
 results[[i]] <- sig_genes %>%
    filter(coef == coefs[i]) %>%
    pull(gene)
  
 names(results)[i] <- coefs[i]
  
  
}


upset(fromList(results))



temp <- results$PostExc_control[ !(results$PostExc_control %in%
                             c(results$W3_copd, 
                               results$W3_control, 
                               results$PostExc_copd))]


mm_results$mixed_model %>%
  filter(model == "nf") %>%
  group_by(model, coef) %>%
  mutate(p.adj = p.adjust(p.val)) %>%
  filter(coef %in% c("timePostExc:copdcopd", 
                     "timePostExc:copdcontrol"), 
         p.adj < 0.05) %>%
  dplyr::select(gene, coef, estimate, p.val, p.adj)
  
  
  
  filter(gene %in% c(temp), 
         coef == "timePostExc:copdcopd", 
         model == "nf") %>%
  ggplot(aes(1, -log10(p.val), size = -estimate)) + geom_point()
















