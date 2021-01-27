library(tidyverse)
library(data.table)
library(dtplyr)
library(readxl)

### Prepare sequencing data for count modeling ###

# This script reads raw counts and meta data, prepares a dgeList,
# filter expressions and saves the list. 

# The files are duplicate sequencing runs. Two dge lists are prepared, 
# one with the duplicate retained and one with the duplicates summarized per
# sample.

# Note that the script produces 13360 and 15093 genes after filtering and 
# 560 and 280 samples for the duplictae and summarised approach respectively 
# (see dim(dge)/dim(dge_sum) below). Gene counts are affected by filtering and summation.


# Read count data 
rc <- read_delim("./data/seq_data/copd_counts.csv", delim = ",")


## Calculate how many zero counts there are
zero_sum <- function(x) {
        s = sum(x == 0)
        return(s)
}

zs <- apply(rc[,-1],1, zero_sum)



# Read sample data 
samples <- read_excel("./data/subject_sample.xlsx", na = "NA") %>%
        mutate(time = Timepoint, 
               time= gsub("3W", "W3", time), 
               subject = as.character(FP), 
               ex.nr = as.character(ex.nr)) %>%
                dplyr::select(-Timepoint, -FP) %>%
        print()



files <- colnames(rc)[-1]

# Renames files to make more workable
files <- gsub("PreExc", "_PreExc_", files)
files <- gsub("PostExc", "_PostExc_", files)
files <- gsub("PreSupp", "_PreSupp_", files)
files <- gsub("3W", "_W3_", files)
files <- gsub("VLL", "VLL_", files)
files <- gsub("VLR", "VLR_", files)
files <- gsub("-", "_", files)


# set colnames on the count table
colnames(rc) <- c("gene", files)


# The counts are derived from duplicate sequencing runs. An alternative dge list is prepared to 
# where counts are summarized across duplicates 

rc <- as.data.table(rc)


rc <- melt(rc, id.vars = "gene", 
              measure.vars = 2:561)


rc[, c("misc.nr", 
           "subject", 
           "time", 
           "leg", 
           "ex.nr", 
           "sample", 
           "lane") := tstrsplit(variable, "_", fixed = TRUE)]



 lazy_dt(rc) %>%
         filter(subject %in% c("102", "144", "105")) %>%

         group_by(subject, gene, time, leg, ex.nr) %>%
         mutate(min = min(value), 
                sum = sum(value)) %>%
         group_by(gene) %>%
         summarise(diff = mean(sum) - mean(min), 
                   min = mean(min)) %>%
  
         filter( min == 0 & diff > 5) %>%
         as_tibble()
         
         distinct(gene) %>%
         pull(gene)
         
         

         ggplot(aes(diff, gene)) + geom_point()
         
         
         print()
         
         
         
        
        filter(gene %in% temp.gene[1:3], subject %in% c("102", "103")) %>%
         
      as_tibble() %>%

         
        ggplot(aes(time, value, group = paste(gene, subject, leg), color = gene)) + 
                geom_line() + 
         geom_point(position = position_jitter()) +
         facet_wrap(~ subject)
        



ct_sum <- rc[, .(sum = sum(value), min = min(value)), by = .(gene,
                                                 subject, 
                                                time, 
                                                 leg, 
                                                 ex.nr)]


temp <- ct_sum[gene %in% genes & min == 0 & sum > 10]




temp[gene %in% genes]




lazy_dt(ct_sum) %>%
        filter(sum > 40, min == 0) %>%
        as_tibble() %>%

        ggplot(aes( sum,subject, color = gene)) + geom_point() 




ct_sum <- rc %>%
       separate(variable, into = c("misc.nr", 
                                "subject", 
                                "time", 
                                "leg", 
                                "ex.nr", 
                                "sample", 
                                "lane")) %>%
        print()
        
        
        group_by(gene, subject, time, leg, ex.nr) %>%
        summarise(min = min(count), 
                  count = sum(count)) %>%
        ungroup() %>%
        print()
        






## This takes some time...










ct_sum <- rc %>%
        
        pivot_longer(names_to = "file",
                     values_to = "count", 
                     cols = 2:561) %>%
        
        separate(file, into = c("misc.nr", 
                                 "subject", 
                                 "time", 
                                 "leg", 
                                 "ex.nr", 
                                 "sample", 
                                 "lane")) %>%
        
        group_by(gene, subject, time, leg, ex.nr) %>%
        summarise(min = min(count), 
                  count = sum(count)) %>%
        ungroup() %>%
        print()


ct_matrix <- ct_sum %>%
        mutate(file = paste0("S", subject, "_", time,"_", leg, "_E", ex.nr)) %>%
        dplyr::select(gene, file, count) %>%
        pivot_wider(names_from = file, values_from = count) %>%
        data.frame() %>%
        print()


ct_mat <- as.matrix(ct_matrix[, -1])
rownames(ct_mat) <- ct_matrix[, 1]

# creates a samples data frame for the combined data
design.df2 <- data.frame(files = colnames(ct_mat)) %>%
        separate(files, into = c("subject", 
                                 "time", 
                                 "leg", 
                                 "ex.nr")) %>%
        mutate(subject = gsub("S", "", subject), 
               ex.nr = gsub("E", "", ex.nr)) %>%
        inner_join(samples) %>%
        print()



# Creates a sample data frame for the combined data
design.df <- data.frame(files = files) %>%
        separate(files, into = c("misc.nr", 
                                 "subject", 
                                 "time", 
                                 "leg", 
                                 "ex.nr", 
                                 "sample", 
                                 "lane")) %>%
        inner_join(samples) %>%
        print()


########### Create DEGlist of all count methods #################

###### Filter and normalization ###########


nonzero    <- rowSums(rc[,-1]) > 0
sum_nonzero <- rowSums(ct_mat) > 0


# Remove zero counts
ct     <- rc[nonzero,]
sum_ct <- ct_mat[sum_nonzero, ]


# Create matrix 
y     <- as.matrix(data.frame(ct)[,-1])
# Set rownames
rownames(y)   <- data.frame(ct)[,1]

# Create dge list
dge <- DGEList(y)
dge_sum <- DGEList(sum_ct)



###### Set samples, combine meta data  #######

## Duplicate dge list ##
dge$samples <- dge$samples %>%
        mutate(files = rownames(.), 
               samples = files) %>%
        separate(files, into = c("misc.nr", 
                                 "subject", 
                                 "time", 
                                 "leg", 
                                 "ex.nr", 
                                 "sample", 
                                 "lane")) %>%
        mutate(misc.nr = gsub("X", "", misc.nr)) %>%
        inner_join(design.df, by = c("misc.nr", "subject", "time", "leg", "ex.nr", "sample", "lane")) %>%
        `rownames<-` (.$samples) %>%
        dplyr::select(-samples) %>%
        dplyr::select(group, lib.size, 
                      norm.factors, 
                      subject, 
                      time, 
                      leg, 
                      ex.nr, 
                      sample,
                      lane,
                      RM_leg,
                      wet.weight, 
                      RQI,
                      treat) %>%
        data.frame()
        
## Sum dge list
dge_sum$samples <- dge_sum$samples %>%
        mutate(files = rownames(.), 
               samples = files) %>%
        separate(files, into = c("subject", 
                                 "time", 
                                 "leg", 
                                 "ex.nr")) %>%
        mutate(subject = gsub("S", "", subject), 
               ex.nr = gsub("E", "", ex.nr)) %>%

        inner_join(design.df2, by = c("subject", "time", "leg", "ex.nr")) %>%
        `rownames<-` (.$samples) %>%
        dplyr::select(-samples) %>%
        dplyr::select(group, lib.size, 
                      norm.factors, 
                      subject, 
                      time, 
                      leg, 
                      ex.nr, 
                      RM_leg,
                      wet.weight, 
                      RQI,
                      treat) %>%
        data.frame()





#### Filter away low expression genes #############

# Create logCPM-plot from all methods #



# Calculate L and M values used for cutoff in each method #
# Note: This cutoff is suggested mased on library sizes, but 
# filterByExpr (used below) uses lower cutoff.

LM.fun <- function(lib){
        
        L <- mean(lib) * 1e-6
        M <- median(lib) * 1e-6
        
        return(data.frame(L = L, M = M))
        
}

cutoff <- rbind(data.frame(LM.fun(dge$samples$lib.size)), 
                data.frame(LM.fun(dge_sum$samples$lib.size))) %>%
        mutate(method = c("duplicates", "sum"), 
               cutoff = log2(10/M + 2/L)) %>%
        print()


# Calculate avg CPM
aCPM <- rbind(data.frame(logCPM = aveLogCPM(dge), method = "duplicates"), 
              data.frame(logCPM = aveLogCPM(dge_sum), method = "sum"))


# Plot of distribution of average log2-CPM. The distribution should be bimodal with a low-abundance peak and
# high abundance peak representing expressed genes. The threshold should separate the two peaks

# Density plot prior to filtering  

denPlot1 <- aCPM %>%
        ggplot(aes(x = logCPM)) + 
        geom_segment(data = data.frame(cutoff), aes(x = cutoff, y = 0, xend = cutoff, yend = 1), 
                     color="red", linetype="dashed") +
        geom_density(alpha = 0.2, stat = "density") + scale_x_continuous(limits = c(-5, 15), expand = c(0,0)) +
        facet_wrap(~ method)



# Save non filtered dge list

saveRDS(dge, file = "./data/derivedData/dge_lists/dge_nonfiltered.RDS")
saveRDS(dge_sum, file = "./data/derivedData/dge_lists/dge_sum_nonfiltered.RDS")



###### Filtering by expression edgeR ################
# Note: this also adds the design matrix to the list.
# The design matrix must be explicitly retrieved in downstream functions,
# e.g. design = dge$design


dge$design <- model.matrix(~ time * treat, data = dge$samples)
dge_sum$design <- model.matrix(~ time * treat, data = dge_sum$samples)


keep.exprs <- filterByExpr(dge, design = dge$design)
keep.exprs_sum <- filterByExpr(dge_sum, design = dge_sum$design)

dge <- dge[keep.exprs, ,keep.lib.sizes = FALSE]
dge_sum <- dge_sum[keep.exprs_sum, ,keep.lib.sizes = FALSE]


aCPM.postfilter <- rbind(data.frame(logCPM = aveLogCPM(dge), method = "duplicates"), 
                         data.frame(logCPM = aveLogCPM(dge_sum), method = "sum"))

logCPM_data <- rbind(data.frame(aCPM, filtering = "pre"), 
                     data.frame(aCPM.postfilter, filtering = "post"))


## Plot of pre and post filtered data
logCPM_data %>%
        ggplot(aes(x = logCPM, fill = filtering)) + 
        geom_segment(data = data.frame(cutoff), 
                     aes(x = cutoff, y = 0, xend = cutoff, yend = 1, fill = NULL), 
                     color="red", linetype="dashed") +
        geom_density(alpha = 0.2, stat = "density") +
        scale_x_continuous(limits = c(-5, 15), expand = c(0,0)) + 
        facet_wrap(~ method)


######## Calculate normalization factor using the TMM method ################

dge <- calcNormFactors(dge, method = "TMM")
dge_sum <- calcNormFactors(dge_sum, method = "TMM")

########## Save dge list ###############

saveRDS(dge, file = "./data/derivedData/dge_lists/dge.RDS")
saveRDS(dge_sum, file = "./data/derivedData/dge_lists/dge_sum.RDS")


dim(dge$counts)
dim(dge_sum$counts)










