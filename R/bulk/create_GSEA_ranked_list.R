##------------------------------------------
## Author: Matina Fragkogianni
## Date: 12-May-2020
##
## A function that creates a ranked file for GSEA analysis 
##
## @param table ,differential expression table
## @param ranked_table, a ranked table
##
##------------------------------------------
create_GSEA_ranked_list <- function(table){
  
  #create ranks file
  ranks_RNAseq = sign(table$logFC) * -log10(table$PValue)
  ranks_RNAseq <- cbind(table$allGenes, ranks_RNAseq)
  colnames(ranks_RNAseq) <- c("GeneName","rank")
  
  #sort ranks in decreasing order
  ranks_RNAseq <- ranks_RNAseq[order(as.numeric(ranks_RNAseq[,2]),decreasing = TRUE),]
  return(ranks_RNAseq)
}