---
title: "Comparable biological and functional adaptations for COPD and healthy participants after 13 weeks of one-legged resistance training"
output:
  html_document:
    keep_md: yes
editor_options: 
  chunk_output_type: console
bibliography: resources/bibliography.bib
---


Correspondence:
Knut Sindre Mølmen,
Section for Health and Exercise Physiology,
Inland Norway University of Applied Sciences,
P.O. Box. 422, 2604 Lillehammer.
E-mail: knut.sindre.molmen@inn.no

Maintainer of this repository: 
Daniel Hammarström, 
Section for Health and Exercise Physiology,
Inland Norway University of Applied Sciences. 
E-mail: daniel.hammarstrom@inn.no
</font>

## Repository overview

Data (raw counts and derived data) are stored in the data folder together with participant characteristics. All R-code can be found in the R-folder. Figures (code and output) are stored in the figures folder.

## Description of R-files

The files are described in the order of their execution.

`seq_prepare.R` Prepares raw count tables for modeling. After basic filtering, dge-lists are stored and used in subsequent scripts.

`mixed_model.R` Models counts for analyses, returns data frames with coefficients, estimates and p-values together with model diagnostics (uniformity tests of simulated residuals).

`cerno_analysis.R` Performs rank-based enrichment tests (non-directional) of gene sets (gene ontology) for coefficients of interest. 

`gene_ontology.R` Similar to cerno, this script performs enrichment tests but a directional gased on fgsea. 

`gene_set_descriptives.R` Calculates descriptive statistic for each gene ontology gene set in (e.g. number of genes with nominal *P*-values < 0.05 etc.)

`kegg_hallmark_gsea_cerno.R` Performs gsea and cerno tests using alternative gene sets (hallmark and KEGG).

`ora_supplementation.R` Performs over representation analysis of differentially expressed genes after supplementation period between tretment and placebo groups.

`lib_fun.R` executed in each script, contains all packages and some functions used.

## Methods 






Gene counts were modeled using negative binomial GLMM with the total library size modeled as a fixed effect [@RN2366] together with sex and COPD status (COPD and healthy controls).
For analyses of the effect of COPD status over time, differential expression was evaluated using GLMMs containing the interaction between time and COPD status. In all models, a single random effect was used, giving each participant an individual intercept. Models were iteratively fitted using glmmTMB [@RN2626].
Model adequacy was assessed for each model fit by assessing uniformity of simulated residuals [@DHARMa].
A total of 15093 were included in the RNA-seq data set after initial filtering,
5.1% were removed due to violation of the uniformity assumption (*P*-value < 0.05).
Genes were identified as differentially expressed when the absolute log2 fold-change was greater than 0.5 and the adjusted P-value (false discovery rate adjusted per model coefficient) was below 5 %. Enrichment analyses of gene sets were performed using two approaches. First, a non-parametric rank test [@RN2439; @RN2438]
was performed based on gene-specific minimum significant differences (MSD).
MSD was defined as the lower limit of the 95 % confidence interval (CI, based on estimated standard errors) around the log fold-change (FC) when log(FC) > 0 and the negative inverse of the upper 95 % CI when log(FC) < 0. Genes with MSD < 0 were further ranked based on P-values. The rank test assessed non-directional changes in gene sets.
Second, gene set enrichment analysis (GSEA) [@RN2432]
was performed to quantify directional regulation of the gene set. GSEA was performed using the fgsea package [@RN2434]
with -log10(p-values) &#215; log2(fold-change) acting as the gene level metric [@RN2627].
Consensus results were given higher importance. Gene ontology (biological process, cellular component and molecular function), as well as Hallmark and KEGG gene sets were retrieved from the molecular signature database (version 7.1) [@RN2436].
Overview of enrichment analyses with exact p-values are presented in supplementary tables. 

## Figure caption

**Figure 1. Transcriptome analyses.** (a), Differences in gene expression between COPD and Healthy at baseline with leading edge genes highlighted from two gene sets identified with gene enrichment analyses (see Table X). Average fold-differences (COPD / Healthy) of genes contributing to baseline differences in from each gene set are shown as individual data points in **b**, violin plots show the distribution of all leading edge genes from each gene set. The average development of each gene set over time are shown in C, the dotted line indicates the mean fold-change of all genes contributing to differences in change over time between COPD status (**c**).         


## Refereces



