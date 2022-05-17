library(tidyverse)
library(edgeR)
library(reshape2)
library(matrixStats)
library(pheatmap)
library(wesanderson)
library(org.Mm.eg.db) 
library(stats)
library(ggplot2)
library(ggrepel)
library(ggthemr)

source(file = "/juno/work/geissmann/data/R_projects/Adipose_tissue_macs/Source/remove.dots.R")

load("/juno/work/geissmann/data/R_projects/bulk_microglia/Results/featureCounts/bulk_microglia.matrix.RData")


# ##---------------------------------------------
# # Load data and make Phenotable
# ##---------------------------------------------

all_counts <- fc$counts
colnames(all_counts)

counts_nuc <- all_counts[,41:48]

list <- lapply(X = colnames(counts_nuc), FUN = strsplit, split = "[.]" )
list <- lapply(list, unlist)

# rename colnames
for(i in 1:length(list)){
    colnames(counts_nuc)[i] <- paste(list[[i]][11], list[[i]][12], list[[i]][13], sep = '_') #sometimes 10
}

# make pheno table
pheno <- data.frame(sample_names = colnames(counts_nuc))
for(i in 1:length(list)){
    pheno$condition[i] <- list[[i]][12]
}


for(i in 1:length(list)){
    print(list[[i]][13])
    list2 <- lapply(X = list[[i]][13], FUN = strsplit, split = "[_]" )
    list2 <- lapply(list2, unlist)
    pheno$B_region[i] <- list2[[1]][1]
}

#pheno$age <- ifelse(grepl(pattern = "2mo", x = pheno$sample_names), "2mo", "ES")
pheno$OG_name <- colnames(all_counts[,41:48])


# ##---------------------------------------------
# # Exploratory analysis for the extra samples
# ##---------------------------------------------

# Create a dataframe with the cpm tranformed values to be plotted
cor_table <- c()
for(i in seq(1, dim(pheno)[1], by=2)){
    df <- data.frame(replicate1 = cpm(counts_nuc[,pheno$sample_names[i]],log=TRUE,prior.count=1),
                     replicate2 = cpm(counts_nuc[,pheno$sample_names[i+1]],log=TRUE,prior.count=1))
    
    # Correlation testing with p-values, using Spearman correlation because data are not normally distributed
    rho = cor.test(x = df$replicate1, y = df$replicate2, method = "spearman", exact = F)
    cor_table <- append(cor_table, rho$estimate)
    
}

##---------------------------------------------
# Merge technical replicates
##---------------------------------------------
counts_nuc[,"M2617_VE_BS_IGO_11769_24_S64_L002_R1_001"] = as.matrix(counts_nuc[,"M2617_VE_BS_IGO_11769_24_S64_L002_R1_001"] + counts_nuc[,"M2617_VE_BS_IGO_11769_24_S64_L003_R1_001"])
counts_nuc[,"M2617_VE_Cort_IGO_11769_23_S63_L002_R1_001"] = as.matrix(counts_nuc[,"M2617_VE_Cort_IGO_11769_23_S63_L002_R1_001"] + counts_nuc[,"M2617_VE_Cort_IGO_11769_23_S63_L003_R1_001"])
counts_nuc[,"M2618_Ctrl_BS_IGO_11769_22_S62_L002_R1_001"] = as.matrix(counts_nuc[,"M2618_Ctrl_BS_IGO_11769_22_S62_L002_R1_001"] + counts_nuc[,"M2618_Ctrl_BS_IGO_11769_22_S62_L003_R1_001"])
counts_nuc[,"M2618_Ctrl_Cort_IGO_11769_21_S61_L002_R1_001"] = as.matrix(counts_nuc[,"M2618_Ctrl_Cort_IGO_11769_21_S61_L002_R1_001"] + counts_nuc[,"M2618_Ctrl_Cort_IGO_11769_21_S61_L003_R1_001"])



# remove the technical replicates from the count matrix
counts_nuc <- counts_nuc[,-which(colnames(counts_nuc) %in% c('M2617_VE_BS_IGO_11769_24_S64_L003_R1_001',
                                                 'M2617_VE_Cort_IGO_11769_23_S63_L003_R1_001',
                                                 'M2618_Ctrl_BS_IGO_11769_22_S62_L003_R1_001',
                                                 'M2618_Ctrl_Cort_IGO_11769_21_S61_L003_R1_001'))]

colnames(counts_nuc)
dim(counts_nuc)

# remove the technical replicates from the phenotable
pheno <- pheno[-which(pheno$sample_names %in% c('M2617_VE_BS_IGO_11769_24_S64_L003_R1_001',
                                                'M2617_VE_Cort_IGO_11769_23_S63_L003_R1_001',
                                                'M2618_Ctrl_BS_IGO_11769_22_S62_L003_R1_001',
                                                'M2618_Ctrl_Cort_IGO_11769_21_S61_L003_R1_001')),]


##---------------------------------------------
# Data importing and normalization
##---------------------------------------------
#make a DGE list
d <- DGEList(counts = counts_nuc, samples = pheno, group = pheno$condition)
d$samples$group <- relevel(x = d$samples$group, ref = "Ctrl")
d$samples$B_region <- as.factor(d$samples$B_region)
d$samples$conRegion <- as.factor(paste(d$samples$group, d$samples$B_region, sep = '.'))
#d$samples$conRegionAge <- as.factor(paste(d$samples$group, d$samples$B_region, d$samples$age, sep = '.'))
d$samples$Age <- as.factor('2mo')


# filtering genes that have enough counts
keep <- filterByExpr(d)
table(keep)
d <- d[keep, , keep.lib.sizes=FALSE]


#TMM normalization
d <- calcNormFactors(d, method = "TMM")

##------------------------------
## Library size plot
##------------------------------
ggthemr('fresh')
g <- ggplot(data = d$samples, aes(x = sample_names, y = lib.size, fill = group)) + geom_bar(stat="identity") 
g + theme(axis.text.x = element_text(angle=45, hjust=1)) 
g

##------------------------------------------------------------
## Log transformation and barplot plot of normalized counts
##------------------------------------------------------------
rownames(d$counts) <- remove.dots(d$counts)
lcpm <- cpm(d, log = TRUE, prior.count = 1)
data = data.frame(t(lcpm), group = d$samples$group, sample_names = d$samples$sample_names)
melted = melt(data = data, id.vars = c("group", "sample_names"))
ggplot(data = melted, aes(x = sample_names, y = value, fill = group)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ylab("Log2 normalised counts per million")

plotDensities(lcpm, legend = FALSE, main = "After filtering")
abline(v = 0, lty = 3)


##-----------------------------------
## Correlation heatmap of samples
##-----------------------------------
df <- data.frame(B_region = d$samples$B_region, Group = d$samples$group)
rownames(df) <- colnames(lcpm)
names(df) <- c("B_region", 'Group') 

Group        <-  wesanderson::wes_palette(name = "Darjeeling1", n = 2)
names(Group) <- levels(df$Group)

B_region        <- wesanderson::wes_palette(name = "Darjeeling2", n = 2)
names(B_region) <- levels(df$B_region)
anno_colors <- list(B_region = B_region, Group = Group)

cor <- cor(lcpm, method = "pearson")
pheatmap(cor, annotation_col = df, annotation_colors = anno_colors,clustering_method = "average", show_colnames = F, show_rownames = F)

##--------------------
## PCA plot analysis
##--------------------

d$samples

df_pca <- prcomp(t(lcpm))
print(summary(df_pca))
df_out <- as.data.frame(df_pca$x)
df_out$sampleConditions <- d$samples$group
df_out$age <- d$samples$age
df_out$B_region <- d$samples$B_region


theme <- theme(panel.background = element_blank(),panel.border = element_rect(fill = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               strip.background=element_blank(),axis.text.x=element_text(colour="black"),
               axis.text.y = element_text(colour="black"), axis.text = element_text(size = 12), axis.ticks = element_line(colour = "black"),
               plot.margin = unit(c(1,1,1,1),"line"), legend.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 12, face = "bold")
               ,axis.title.y = element_text(color = "black", size = 12, face = "bold")) 
theme_update(plot.title = element_text(hjust = 0.5, face = "bold"))


p <- ggplot(df_out,aes(x = PC1, y = PC2, color = sampleConditions, shape = B_region))
p <- p + geom_point(size = 3, stroke = 0.8)+ theme + scale_color_manual(values =  c('#1339F4', '#C70039')) #+ xlab(percentage[1]) + ylab(percentage[2]) 
p 


##--------------------------------------
## Differential expression analysis 
##--------------------------------------

group <- factor(d$samples$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

#design matrix
#combination_B_region <- factor(d$samples$conPop)
#design <- model.matrix(~0 + combination_B_region)
#colnames(design) <- levels(combination_B_region)

# estimate dispersion
d <- estimateDisp(d, design, robust = TRUE)
d$common.dispersion
plotBCV(d, main="BCV Plot")

cont.matrix <- makeContrasts(
    #CESvsVEES = VE.ES-Ctrl.ES,
    #C2mvsVE2m = VE.2mo-Ctrl.2mo,
    CvsVE = VE-Ctrl,
    #CvsVE_BS = VE.BS-Ctrl.BS,
    #CvsVE_Cort = VE.Cort-Ctrl.Cort,
    levels = design)


fit <- glmQLFit(d, design, robust = T)
qlf <- glmQLFTest(fit, contrast=cont.matrix[,"CvsVE"])
tt <- topTags(qlf, n = 99999999)
table <- tt$table
summary(decideTests(qlf, lfc = 1.5, p.value = 0.05))


#### Filter the gene table to identify significant genes
sigGenes <- tt[table$FDR <= 0.05,] #get seperataly the significant genes
SigGenesUpDown <- tt[table$FDR <= 0.05 & (table$logFC >= 1.5  | table$logFC <= -1.5),] #get seperataly the significant genes
UpGenes <- sigGenes[sigGenes$table$logFC >= 1.5, ]
DownGenes <- sigGenes[sigGenes$table$logFC <= -1.5, ]

dim(SigGenesUpDown)
dim(UpGenes)
dim(DownGenes)

allGenes <- mapIds(org.Mm.eg.db, keys = rownames(table), keytype = "ENSEMBL", column="SYMBOL")
allGenes <- as.data.frame(allGenes)
testVlookup <- merge(table, allGenes, by.x="row.names", by.y="row.names") 



##-------------------------------------------------
## Heatmap of common genes between Human and BS 2mo 
##-------------------------------------------------
UP_human_mouse_common_2mo <- read.delim("/juno/work/geissmann/data/R_projects/bulk_microglia/Results/Human_Mouse_Common_Genes/UP_human_mouse_common_2mo.txt")
matrix <- lcpm[which(rownames(lcpm) %in% UP_human_mouse_common_2mo$Mouse.Gene_ID),] #sig_up$allGene
#matrix <- lcpm[which(toupper(rownames(lcpm)) %in% toupper(genes_up)),] #genes_up
dim(matrix)

df <- data.frame(Age = d$samples$Age,B_region = d$samples$B_region, Group = d$samples$group)
rownames(df) <- colnames(lcpm)
names(df) <- c('Age',"B_region", 'Group') 

Group        <-  c('#1339F4', '#C70039')
names(Group) <- levels(df$Group)

B_region        <- c('#22BE6E', '#FFC300')
names(B_region) <- levels(df$B_region)

Age       <- c('#282409')
names(Age) <- levels(as.factor(df$Age))

anno_colors <- list(Age = Age, B_region = B_region, Group = Group)

phenotable_ordered <- pheno %>% arrange(condition, B_region)
matrix_ordered <- matrix[,phenotable_ordered$sample_names]


#Gene names
gene_names <- mapIds(org.Mm.eg.db, keys = rownames(matrix_ordered), keytype = "ENSEMBL", column="SYMBOL")
gene_names <- as.data.frame(gene_names)
rownames(matrix_ordered) <- gene_names$gene_names

pheatmap(matrix_ordered, annotation_col = df, annotation_colors = anno_colors,  scale = "row", 
         show_rownames = T, 
         show_colnames = F, cluster_cols = F, cluster_rows = F)

##-----------------------------------
## Neuroinflammation analysis
##-----------------------------------

genes_up <- c('CHIT1',
              'IL1B',
              'CHI3L1',
              'CXCL10',
              'LYZ',
              'NLRC4',
              'CXCL11',
              'PLA2G4A',
              'SCIN',
              'NCF1',
              'TREM2',
              'C1QC',
              'C1QB',
              'CFI',
              'IL1A',
              'CX3CR1',
              'CSF1R',
              'SOCS3',
              'CTSS',
              'TYROBP',
              'MMP9',
              'MSR1',
              'CTSC',
              'IFITM1')

##-----------------------------------
## Heatmap of Neuroinflammatory genes 
##-----------------------------------
#Gene names
gene_names <- mapIds(org.Mm.eg.db, keys = rownames(lcpm), keytype = "ENSEMBL", column="SYMBOL")
gene_names <- as.data.frame(gene_names)
#testVlookup <- merge(matrix, gene_names, by.x="row.names", by.y="row.names") 
rownames(lcpm) <- gene_names$gene_names

#matrix <- lcpm[which(toupper(rownames(lcpm)) %in% toupper(sig_up$allGenes),] #sig_up$allGene
matrix <- lcpm[which(toupper(rownames(lcpm)) %in% toupper(genes_up)),] #genes_up
dim(matrix)

df <- data.frame(B_region = d$samples$B_region, Group = d$samples$group)
rownames(df) <- colnames(lcpm)
names(df) <- c("B_region", 'Group') 

Group        <-  c('#1339F4', '#C70039')
names(Group) <- levels(df$Group)

B_region        <- c('#22BE6E', '#FFC300')
names(B_region) <- levels(df$B_region)
anno_colors <- list(B_region = B_region, Group = Group)

phenotable_ordered <- pheno %>% arrange(condition, B_region)
matrix_ordered <- matrix[,phenotable_ordered$sample_names]

gene_order <- c('Cxcl10', 'Il1b', 'Msr1', 'Socs3', 'Ctss', 'Pla2g4a','Il1a', 'C1qb', 'Ncf1','Tyrobp', 'Cx3cr1', 'Trem2', 'Csf1r', "Mmp9", "Ctsc", "Ifitm1", "Nlrc4")
matrix_ordered <- matrix_ordered[gene_order,]

pheatmap(matrix_ordered, annotation_col = df, annotation_colors = anno_colors,  scale = "row", 
         show_rownames = T, 
         show_colnames = F, cluster_cols = F, cluster_rows = F)



