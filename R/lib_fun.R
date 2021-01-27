# Libraries for analysis #

# BiocManager::install("biomaRt")


library(piano)



## Bioconductor packages
library(org.Hs.eg.db)
library(msigdbr)

library(fgsea)

library(edgeR)
library(limma)


library(qvalue)
library(preprocessCore)

library(CePa)

library(topGO)
library(GOSim)

library(GOSemSim)
library(clusterProfiler)
library(DOSE)

library(DHARMa)

library(biomaRt)


# Cran packages
library(tmod)
library(RColorBrewer)
# library(UpSetR)
library(ggtext)


library(ggraph)
# library(enrichplot)
library(tidygraph)

library(tidyverse)
library(ggrepel)

library(readxl)
library(nlme)
library(emmeans)
# library(feather)

library(glmmTMB)
library(lme4)
# library(mgcv)

library(svMisc)
library(foreach)

library(doSNOW)
library(parallel)



### For plotting 
library(cowplot)



############ Functions ##############



# Makes first letter of string to capital 

####### Makes first letter of string to capital ###################


firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}



##### Function for splitting up strings with newline #####

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



