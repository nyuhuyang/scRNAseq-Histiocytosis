##------------------------------------------
## Author: Matina Fragkogianni
## Date: 115-09-2021
##
## A function to read in featureCounts results and create a counts table and a phenodata matrix 
##
## @param fc is the Rdata object
## @return a list containing counts and phenodata
##
##------------------------------------------

import.featureCount.data <- function(fc, position){
    
    counts <- fc$counts
    colnames(counts)
    
    list <- lapply(X = colnames(counts), FUN = strsplit, split = "[.]" )
    list <- lapply(list, unlist)
    
    # rename colnames
    for(i in 1:length(list)){
        colnames(counts)[i] <- list[[i]][position] #sometimes 10
    }
    
    # make pheno table
    pheno <- data.frame(sample_names = colnames(counts))
    pheno$population <- ifelse(grepl(pattern = "Alpha", x = pheno$sample_names, ignore.case = T), "Alpha", ifelse(grepl(pattern = "delta", x = pheno$sample_names, ignore.case = T), "Delta", ifelse(grepl(pattern = "gamma", x = pheno$sample_names, ignore.case = T), "Gamma", ifelse(grepl(pattern = "Mam", x = pheno$sample_names, ignore.case = T), "Mam","??"))))
    pheno$condition <- ifelse(grepl(pattern = "py", x = pheno$sample_names), "PyMT", "WT")
    pheno$conPop <- paste(pheno$population, pheno$condition, sep = ".")
    
    return(list(counts = counts, pheno = pheno))
}