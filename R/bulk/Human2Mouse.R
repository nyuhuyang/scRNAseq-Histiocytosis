Human2Mouse <- function(HumanGenes){
    
    human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    genesHumantoMouse <- biomaRt::getLDS(attributes = c("ensembl_gene_id","hgnc_symbol"), 
                                        filters = "ensembl_gene_id", 
                                        values = HumanGenes, 
                                        mart = human, 
                                        attributesL = c("ensembl_gene_id", "mgi_symbol"), 
                                        martL = mouse, 
                                        uniqueRows = TRUE)
    
    colnames(genesHumantoMouse) <- c("Human.Gene_ID", "HGNC","Mouse.Gene_ID", "MGI")
    
    return(genesHumantoMouse) 
    
}