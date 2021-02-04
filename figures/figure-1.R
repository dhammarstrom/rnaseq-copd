##-------------------------------------
## figure-1.R
##
## Title: Figure 1. Baseline differences between COPD and CONTROLS
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

# Librarises
source("./R/lib_fun.R")





# Themes, colors, functions etc. #######################





# Color scale 

col.scale <- c("#a6cee3",
                "#1f78b4",
                "#b2df8a",
                "#33a02c")




# Plotting theme 

plot_theme <- function(){
  
  theme_bw() +
    theme(axis.text = element_text(size = 8), 
          axis.title = element_text(size = 8), 
          legend.title = element_text(size = 8), 
          legend.text = element_text(size = 8), 
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "gray95"))
  
}


new_line_fun <- function(x){
  
  l <- sapply(strsplit(x, " "), length)
  
  if(l > 2){
    
    paste(paste(strsplit(x, " ")[[1]][1:2], collapse = " "), 
          "\n", 
          paste(strsplit(x, " ")[[1]][3:l], collapse = " "), 
          collaspe = "")
    
  } else {
    x
  }
  
}



## Read data for gene expression data individual genes ##########################


mixed_dat <- readRDS(file = "./data/derivedData/mixed_model/mixed_full_copd.RDS")



dat <- mixed_dat %>%
  filter(model == "nf", 
         uniform.p > 0.05, 
         coef == "copdcopd") %>% 
  mutate(p.adj = p.adjust(p.val, method = "fdr"), 
         log2fc = estimate/log(2)) %>%
  print()


# named vectors of gene-level statistics
gene_name_convert <- clusterProfiler::bitr(unique(dat$gene), fromType = "ENSEMBL",
                                           toType = c( "ENTREZID", "SYMBOL"),
                                           OrgDb = org.Hs.eg.db)

# df with symbols
symbols <- dat %>%
  left_join(gene_name_convert %>%
               mutate(gene  = ENSEMBL) %>%
               dplyr::select(gene, SYMBOL)) 




########################## Genes category plot ####################


# Load gene set descriptives
gene_sets <- readRDS(file = "./data/derivedData/geneset_descriptives.RDS")

unique(gene_sets$coefs)

gene_sets %>%
  filter(ID %in% c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", 
                   "HALLMARK_MYOGENESIS"), 
         coefs %in% c("nf_copdcopd", 
                      "nf_timeW3:copdcopd", 
                      "nf_timePostExc:copdcopd"))


# Load leading edge genes from gsea analysis

full_gsea_df <- readRDS(file = "./data/derivedData/fgsea_full.RDS")

gsea <- full_gsea_df %>%
  dplyr::select(ID = pathway,
                leadingEdge,
                coef) %>%
  filter(ID %in% c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", 
                   "HALLMARK_MYOGENESIS"), 
         coef %in% c("nf_copdcopd", 
                      "nf_timeW3:copdcopd", 
                      "nf_timePostExc:copdcopd"))


### Extract genes from leading edges

myogenesis_goi_change <- gsea %>%
  filter(ID %in% c("HALLMARK_MYOGENESIS"), 
         coef %in% c("nf_timeW3:copdcopd", 
                     "nf_timePostExc:copdcopd")) %>%
  pull(leadingEdge) %>%
  unlist()

myogenesis_goi_baseline <- gsea %>%
  filter(ID %in% c("HALLMARK_MYOGENESIS"), 
         coef %in% c("nf_copdcopd")) %>%
  pull(leadingEdge) %>%
  unlist()



oxphos_goi_change <- gsea %>%
  filter(ID %in% c("HALLMARK_OXIDATIVE_PHOSPHORYLATION"), 
         coef %in% c("nf_timeW3:copdcopd", 
                     "nf_timePostExc:copdcopd")) %>%
  pull(leadingEdge) %>%
  unlist()

oxphos_goi_baseline <- gsea %>%
  filter(ID %in% c("HALLMARK_OXIDATIVE_PHOSPHORYLATION"), 
         coef %in% c("nf_copdcopd")) %>%
  pull(leadingEdge) %>%
  unlist()




## Volcano plot baseline ##############################


volcano_fig_dat  <- symbols %>%
  
  
  
  #  filter(log2fc < 10, log2fc > -10) %>%
  mutate(gene_set = if_else(SYMBOL %in% myogenesis_goi_baseline, 
                            "myogenesis", 
                            if_else(SYMBOL %in% oxphos_goi_baseline, 
                                    "oxphos", "other")), 
         gene_set = factor(gene_set, levels = c("oxphos", 
                                                "myogenesis", "other"), 
                           labels = c("Oxidative\nphosphorylation", 
                                      "Myogenesis", "other")), 
         sig = if_else(p.adj < 0.05 & abs(log2fc) > 0.5, "sig", "non.sig"), 
         sig2 = if_else(p.adj < 0.02 & abs(log2fc) > 0.55, "sig", "non.sig"), 
         sig.metric = -log10(p.adj) * abs(log2fc)) %>%
  
  print()
  
volcano_fig <- volcano_fig_dat %>%
  filter(gene_set == "other") %>%
  ggplot(aes(log2fc,  -log10(p.adj))) +
  
  # A line to specify p < 0.05
  geom_hline(yintercept = -log10(0.05), 
             lty = 2, 
             alpha = 0.2) +
 
   annotate("richtext", x = -4.8, y = -log10(0.05) + 0.1, 
           label = "*FDR* = 0.05", 
           hjust = 0,
           alpha = 0.5,
           size = 3, 
           fill = "gray95", 
           label.color = NA, # remove background and outline
           label.padding = grid::unit(rep(0, 4), "pt")) +
  
  
  
  geom_point(shape = 21, 
             fill = col.scale[4], 
             alpha = 0.2, 
             size = 2.5) +
  
  geom_point(data = filter(volcano_fig_dat, gene_set != "other"), 
             aes(fill = gene_set), 
             shape = 21, 
             size = 3) +
  
  scale_fill_manual(values = c(col.scale[3], col.scale[1])) +
  
  # Labels for each gene
#   geom_text_repel(data = . %>%
#                     filter(sig.metric > 3.2, 
#                            
#                            log2fc < 0), 
#                   aes(label = SYMBOL), 
#                   size = 2, 
#                   color = "black",
#                   direction = "both", 
#                   nudge_x = -2, 
#                   nudge_y = 1,
#                   alpha = 1, 
#                   segment.alpha = 0.3) +
#   geom_text_repel(data = . %>%
#                     filter(sig.metric > 3.2, 
#                            log2fc > 0), 
#                   aes(label = SYMBOL), 
#                   size = 2, 
#                   color = "black",
#                   direction = "both", 
#                   nudge_x = 2, 
#                   nudge_y = 1,
#                   alpha = 1, 
#                   segment.alpha = 0.3)  +
  
  scale_x_continuous(limits = c(-5,  5), 
                     expand = c(0, 0), 
                     breaks = c(-5, -2.5, 0, 2.5, 5)) +
  scale_y_continuous(limits = c(0, 5), 
                     expand = c(0,0), 
                     breaks = c(0, 2.5, 5)) +
  
  
  
  guides(fill = guide_legend()) +
  
  labs(y = "-Log<sub>10</sub>(adjusted *P*-value)", 
       x = "Log<sub>2</sub> fold-difference (COPD - Healthy)") +
  
  plot_theme() + 
  theme(legend.position = c(0.85, 0.87), 
        legend.margin = margin(0, 0, 0, 0, unit='cm'),
        legend.title = element_blank(), 
        axis.title.x = element_markdown(size = 8), 
        axis.title.y = element_markdown(size = 8))














# Load gsea and cerno statistics

# full gsea
full_gsea_df <- readRDS(file = "./data/derivedData/fgsea_full.RDS")

gsea <- full_gsea_df %>%
  dplyr::select(ID = pathway,
                gsea.padj = padj,
                gsea.p = pval, 
                gsea.es = ES, 
                gsea.nes = NES, 
                coef) 



# full cerno
cerno.test <- readRDS("./data/derivedData/cerno.RDS")

cerno <- cerno.test %>%
  dplyr::select(ID = Title,
                cerno, 
                N1, 
                cerno.p = P.Value,
                cerno.padj = adj.P.Val, 
                go.cat, 
                coef) 



#### Read gene set descriptives  #####
go.descr <- readRDS(file = "./data/derivedData/geneset_descriptives.RDS") %>%
  mutate(coefs = gsub("tissue_offset", "tissue", coefs)) %>%
  separate(coefs, into = c("model", "coef"), sep = "_") %>%
  mutate(coef = gsub("tissue", "tissue_offset", coef)) %>%
  
  mutate(msd.pos = paste0(round((msd.pos/n.total) * 100,1),"%"), 
         msd.pos_inLE = round((msd.pos_inLE/in_LE) * 100, 1), 
         LE = paste0(in_LE, " (", msd.pos_inLE, "%)"), 
         fc = paste0(round(fc.LE, 2), 
                     " [", 
                     round(fc.LE.min,2), 
                     ", ",
                     round(fc.LE.max, 2), 
                     "]")) %>%
  dplyr::select(model, coef, ID, n.total, set.size, msd.pos, LE, fc) %>%
  mutate(set.size = paste0(n.total, " (", set.size, ")")) 



# Top fgsea gene categories
gene_set_stats <- gsea %>%
  inner_join(cerno) %>%

  # Change coef into two variables 
  mutate(coef = gsub("tissue_offset", "tissue", coef)) %>%
  separate(coef, into = c("model", "coef"), sep = "_") %>%
  
  filter(model == "nf", 
         coef %in% c("copdcopd", "timeW3:copdcopd", "timePostExc:copdcopd"), 
         ID %in% c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", 
                   "HALLMARK_MYOGENESIS")) %>%
  print()
  



## Load DE model data

mixed_dat <- readRDS(file = "./data/derivedData/mixed_model/mixed_full_copd.RDS")

## Convert gene ids to symbol

gene_name_convert <- clusterProfiler::bitr(unique(mixed_dat$gene), fromType = "ENSEMBL",
                                           toType = c( "ENTREZID", "SYMBOL"),
                                           OrgDb = org.Hs.eg.db) %>%
  dplyr::select(gene = ENSEMBL, SYMBOL)






## Baseline differences plot #####################


baseline_dat <- mixed_dat %>%
  inner_join(gene_name_convert) %>%
  filter(model == "nf", 
         coef == "copdcopd") %>%
  mutate(gene_set_change = if_else(SYMBOL %in% myogenesis_goi_change, "myogenesis", 
                            if_else(SYMBOL %in% oxphos_goi_change, "oxphos",
                                    "other")),
         
         gene_set_baseline = if_else(SYMBOL %in% myogenesis_goi_baseline, "myogenesis", 
                                     if_else(SYMBOL %in% oxphos_goi_baseline, "oxphos",
                                             "other")), 
         
         gene_set = if_else(SYMBOL %in% c(myogenesis_goi_baseline,
                                          myogenesis_goi_change), "myogenesis", 
                                     if_else(SYMBOL %in% c(oxphos_goi_baseline,
                                                           oxphos_goi_change), "oxphos",
                                             "other"))) %>%
  filter(gene_set != "other") %>%
  mutate(fc = exp(estimate), 
         gene_set = factor(gene_set, levels = c("oxphos", 
                                                "myogenesis"), 
                           labels = c("Oxidative phosphorylation", 
                                      "Myogenesis"))) %>%
  print()
  

# Fold differences at baseline COPD - CONTROL
#

baseline_comp <- baseline_dat %>%
  ggplot(aes(gene_set, fc, fill = gene_set)) + 
  
  geom_hline(yintercept = 1, lty = 2, alpha = 0.2) +
  
  geom_violin(alpha = 0.2) + 
  geom_point(data = filter(baseline_dat, gene_set_baseline != "other"), 
             aes(fill = gene_set), 
             position = position_jitter(width = 0.07), 
             shape = 21, 
             size = 2, 
             alpha = 0.6) +
  
  scale_fill_manual(values = c(col.scale[3], col.scale[2])) +
  
  scale_y_continuous(limits = c(0, 3), 
                     breaks = c(1, 2, 3), 
                     expand = c(0,0)) +
  
  labs(y = "Fold-difference at baseline\n(COPD - Healthy)") +
  
  plot_theme() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "none")





# Create a data frame with all 1 at pre
means <- mixed_dat %>%
  
  filter(model == "nf", 
         interaction == TRUE, 
         uniform.p > 0.05) %>%
  dplyr::select(gene, coef, estimate) %>%
  pivot_wider(names_from = coef, 
              values_from = estimate) %>%
  # Marginal means 
  mutate(control_pre = `(Intercept)`,
         control_w3 = `(Intercept)` + timeW3,
         control_post = `(Intercept)` + timePostExc, 
         copd_pre = control_pre + copdcopd, 
         copd_w3 = control_w3 + copdcopd + `timeW3:copdcopd`, 
         copd_post = control_post + copdcopd + `timePostExc:copdcopd`) %>%
  dplyr::select(gene, control_pre:copd_post) %>%
  
  mutate(copd_w3_fc = exp(copd_w3 - copd_pre), 
         copd_post_fc = exp(copd_post - copd_pre), 
         control_w3_fc = exp(control_w3 - control_pre), 
         control_post_fc = exp(control_post - control_pre)) %>%
  
  
  print()






pre_dat <- means %>%
  inner_join(gene_name_convert) %>%
  mutate(gene_set_change = if_else(SYMBOL %in% myogenesis_goi_change, "myogenesis", 
                                   if_else(SYMBOL %in% oxphos_goi_change, "oxphos",
                                           "other")),
         
         gene_set_baseline = if_else(SYMBOL %in% myogenesis_goi_baseline, "myogenesis", 
                                     if_else(SYMBOL %in% oxphos_goi_baseline, "oxphos",
                                             "other")), 
         
         gene_set = if_else(SYMBOL %in% c(myogenesis_goi_baseline,
                                          myogenesis_goi_change), "myogenesis", 
                            if_else(SYMBOL %in% c(oxphos_goi_baseline,
                                                  oxphos_goi_change), "oxphos",
                                    "other"))) %>%
  filter(gene_set != "other") %>%
  mutate(copd_pre = 1, 
         control_pre = 1) %>%
  dplyr::select(gene, SYMBOL, gene_set,copd_pre, control_pre) %>%
  print()
  

# calculate mean fold change of leading edge genes from change 
# category

cat_means <- means %>%
  inner_join(gene_name_convert) %>%
  dplyr::select(gene, SYMBOL, copd_w3_fc:control_post_fc) %>%
  
  mutate(gene_set_change = if_else(SYMBOL %in% c(myogenesis_goi_change,oxphos_goi_change),
                                   "change", "other")) %>%
  
  inner_join(pre_dat) %>%
  
  pivot_longer(names_to = "group_time", 
               values_to = "fc", 
               cols = c(copd_w3_fc:control_post_fc, copd_pre, control_pre)) %>%
  separate(group_time, into = c("group", "time")) %>%
  filter(gene_set_change == "change") %>%
  
  mutate(time = factor(time, levels = c("pre", 
                                        "w3", 
                                        "post"), 
                       labels = c("Baseline", 
                                  "Week 3&#189;", 
                                  "Post-RT")), 
         group = factor(group, levels = c("copd", "control"), 
                        labels = c("COPD", "Healthy")), 
         gene_set = factor(gene_set, levels = c("oxphos", 
                                                "myogenesis"), 
                           labels = c("Oxidative phosphorylation", 
                                      "Myogenesis"))) %>%
  
  
  group_by(gene_set, time, group) %>%
  summarise(fc = mean(fc)) %>%
  
  print()



fold_change_plot <- means %>%
  inner_join(gene_name_convert) %>%
  dplyr::select(gene, SYMBOL, copd_w3_fc:control_post_fc) %>%

  mutate(gene_set_change = if_else(SYMBOL %in% c(myogenesis_goi_change,oxphos_goi_change),
                                   "change", "other")) %>%

  inner_join(pre_dat) %>%

  pivot_longer(names_to = "group_time", 
               values_to = "fc", 
               cols = c(copd_w3_fc:control_post_fc, copd_pre, control_pre)) %>%
  separate(group_time, into = c("group", "time")) %>%
  
  mutate(time = factor(time, levels = c("pre", 
                                        "w3", 
                                        "post"), 
                       labels = c("Baseline", 
                                  "Week 3&#189;", 
                                  "Post-RT")), 
         group = factor(group, levels = c("copd", "control"), 
                        labels = c("COPD", "Healthy")), 
         gene_set = factor(gene_set, levels = c("oxphos", 
                                                "myogenesis"), 
                           labels = c("Oxidative phosphorylation", 
                                      "Myogenesis"))) %>%
  

  
  ggplot(aes(time, fc, group = paste(SYMBOL, group), 
             color = gene_set, 
             alpha = gene_set_change)) + 
  
  geom_line(size = 0.4) + 
  
  geom_line(data = cat_means, lty = 2, alpha = 1, 
            color = "black",
            aes(group = group)) +
  
  scale_alpha_manual(values = c(.5, .2)) +
  scale_color_manual(values = c(col.scale[3], col.scale[2])) +
  scale_y_continuous(limits = c(0.4, 2.1), 
                     breaks = c(0.5, 1, 1.5, 2), 
                     labels = c("", 1, "", 2)) +

  labs(y = "Fold-change from baseline") +
    
  facet_grid(gene_set ~ group) +
  plot_theme() + 
  theme(axis.text.x  = element_markdown(angle = 45, hjust = 1), 
        axis.title.x = element_blank(), 
        legend.position = "none")
  











###### Save figure 
figure_1 <- plot_grid(
  plot_grid(volcano_fig, baseline_comp,
            align = "hv",
            rel_heights = c(0.5, 0.5), ncol = 1), 
           plot_grid(NULL, fold_change_plot,NULL, ncol = 1, rel_heights = c(0.1, 1, 0.1)),
  ncol = 2, 
                      rel_widths = c(0.55, 0.45)) +
  
  draw_plot_label(label=c("A", "B", "C"),
                  x =   c(0.02, 0.02, 0.57), 
                  y =   c(0.98, 0.48, 0.91),
                  hjust=.5, vjust=.5, size = 14)
                      

# Width of figure = 2x columns 190 mm, 1 column = 90 mm
# height of figure = full page = 23 cm       


ggsave("figures/figure-1.pdf", plot = figure_1, width = 19, height = 24/2, 
       dpi = 600,
       units = "cm", device=cairo_pdf)













