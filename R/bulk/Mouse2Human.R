Mouse2Human <- function(MouseGenes){
  
  
  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesMousetoHuman = biomaRt::getLDS(attributes = c("ensembl_gene_id","mgi_symbol"), 
                             filters = "mgi_symbol", 
                             values = MouseGenes , 
                             mart = mouse, 
                             attributesL = c("ensembl_gene_id", "hgnc_symbol"), 
                             martL = human, 
                             uniqueRows = TRUE)
  
  colnames(genesMousetoHuman) <- c("Mouse.Gene_ID", "MGI", "Human.Gene_ID", "HGNC")
  
  return(genesMousetoHuman) 
  
}
