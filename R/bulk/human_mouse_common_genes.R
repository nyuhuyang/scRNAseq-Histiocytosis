
source(file = "/juno/work/geissmann/data/R_projects/bulk_microglia/source/Human2Mouse.R")

UP_genes_BS_ES <- read.delim("/juno/work/geissmann/data/R_projects/bulk_microglia/Results/DEA_results/ES/UP_genes_BS_ES.txt")
DownGenes_BS_ES <- read.delim("/juno/work/geissmann/data/R_projects/bulk_microglia/Results/DEA_results/ES/DownGenes_BS_ES.txt")

UP_genes_BS_2mo <- read.delim("/juno/work/geissmann/data/R_projects/bulk_microglia/Results/DEA_results/2mo/UP_genes_BS_2mo.txt")
DOWN_genes_BS_2mo <- read.delim("/juno/work/geissmann/data/R_projects/bulk_microglia/Results/DEA_results/2mo/DOWN_genes_BS_2mo.txt")

upregulated_ALL_excludingCERE_gene_names <- read.delim("/juno/work/geissmann/data/R_projects/wholeBrainHuman/Results/DEA/upregulated_ALL_excludingCERE_gene_names_FC_1.5.txt")
dim(upregulated_ALL_excludingCERE_gene_names)
downregulated_ALL_excludingCERE_gene_names <- read.delim("/juno/work/geissmann/data/R_projects/wholeBrainHuman/Results/DEA/downregulated_ALL_excludingCERE_gene_names_FC_1.5.txt")
dim(downregulated_ALL_excludingCERE_gene_names)

########################
#### UP-regulated genes
########################
UP_human2Mouse_ORTHO <- Human2Mouse(HumanGenes = rownames(upregulated_ALL_excludingCERE_gene_names))

#### 2 month old mice &  Human Whole Brain
sprintf("Human orthologs: %d", length(UP_human2Mouse_ORTHO$Mouse.Gene_ID))
sprintf("Mouse: %d", length(UP_genes_BS_2mo$Row.names))
sum(UP_human2Mouse_ORTHO$Mouse.Gene_ID %in% UP_genes_BS_2mo$Row.names)
UP_human_mouse_common_2mo <- UP_human2Mouse_ORTHO[which(UP_human2Mouse_ORTHO$Mouse.Gene_ID %in% UP_genes_BS_2mo$Row.names),]
#write.table(x = UP_human_mouse_common_2mo, file = "Results/Human_Mouse_Common_Genes/UP_human_mouse_common_2mo.txt", sep = "\t")

#### End stage  mice  &  Human Whole Brain
sprintf("Human orthologs: %d", length(UP_human2Mouse_ORTHO$Mouse.Gene_ID))
sprintf("Mouse: %d", length( UP_genes_BS_ES$Row.names))
sum(UP_human2Mouse_ORTHO$Mouse.Gene_ID %in% UP_genes_BS_ES$Row.names)
UP_human_mouse_common_ES <- UP_human2Mouse_ORTHO[which(UP_human2Mouse_ORTHO$Mouse.Gene_ID %in% UP_genes_BS_ES$Row.names),]
#write.table(x = UP_human_mouse_common_ES, file = "Results/Human_Mouse_Common_Genes/UP_human_mouse_common_ES.txt", sep = "\t")


########################
#### Down-regulated genes
########################
DOWN_human2Mouse_ORTHO <- Human2Mouse(HumanGenes = rownames(downregulated_ALL_excludingCERE_gene_names))

#### 2 month old mice  &  Human Whole Brain
sprintf("Human orthologs: %d", length(DOWN_human2Mouse_ORTHO$Mouse.Gene_ID))
sprintf("Mouse: %d", length(DOWN_genes_BS_2mo$Row.names))
sum(DOWN_human2Mouse_ORTHO$Mouse.Gene_ID %in% DOWN_genes_BS_2mo$Row.names)
DOWN_human_mouse_common_2mo <- DOWN_human2Mouse_ORTHO[which(DOWN_human2Mouse_ORTHO$Mouse.Gene_ID %in% DOWN_genes_BS_2mo$Row.names),]

#### End stage  mice  &  Human Whole Brain
sprintf("Human orthologs: %d", length(DOWN_human2Mouse_ORTHO$Mouse.Gene_ID))
sprintf("Mouse: %d", length( DownGenes_BS_ES$Row.names))
sum(DOWN_human2Mouse_ORTHO$Mouse.Gene_ID %in% DownGenes_BS_ES$Row.names)
DOWN_human_mouse_common_ES <- DOWN_human2Mouse_ORTHO[which(DOWN_human2Mouse_ORTHO$Mouse.Gene_ID %in% DownGenes_BS_ES$Row.names),]




