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
library(fgsea)

source(file = "/juno/work/geissmann/data/R_projects/Adipose_tissue_macs/Source/remove.dots.R")

load("/juno/work/geissmann/data/R_projects/bulk_microglia/Results/featureCounts/bulk_microglia.matrix.RData")

# ##---------------------------------------------
# # Load data and make Phenotable
# ##---------------------------------------------

all_counts <- fc$counts
colnames(all_counts)

counts_nuc <- all_counts[,41:48]
counts <- all_counts[,1:40]
colnames(counts)

list <- lapply(X = colnames(counts), FUN = strsplit, split = "[.]" )
list <- lapply(list, unlist)

# rename colnames
for(i in 1:length(list)){
    colnames(counts)[i] <- paste(list[[i]][10], list[[i]][13], sep = '_') #sometimes 10
}

# make pheno table
pheno <- data.frame(sample_names = colnames(counts))
for(i in 1:length(list)){
    pheno$condition[i] <- list[[i]][11] #sometimes 10
}

for(i in 1:length(list)){
    pheno$B_region[i] <- list[[i]][12] #sometimes 10
}

pheno$age <- ifelse(grepl(pattern = "2mo", x = pheno$sample_names), "2mo", "ES")
pheno$OG_name <- colnames(all_counts[,1:40])


# ##---------------------------------------------
# # Exploratory analysis for the extra samples
# ##---------------------------------------------

# Create a dataframe with the cpm tranformed values to be plotted
cor_table <- c()
for(i in seq(1, dim(pheno)[1], by=2)){
    df <- data.frame(replicate1 = cpm(counts[,pheno$sample_names[i]],log=TRUE,prior.count=1),
                 replicate2 = cpm(counts[,pheno$sample_names[i+1]],log=TRUE,prior.count=1))

    # Correlation testing with p-values, using Spearman correlation because data are not normally distributed
    rho = cor.test(x = df$replicate1, y = df$replicate2, method = "spearman", exact = F)
    cor_table <- append(cor_table, rho$estimate)

}

##---------------------------------------------
# Merge technical replicates
##---------------------------------------------
counts[,"M1375_ES_IGO_11769_2_S65_L002_R1_001"] = as.matrix(counts[,"M1375_ES_IGO_11769_2_S65_L002_R1_001"] + counts[,"M1375_ES_IGO_11769_2_S65_L003_R1_001"])
counts[,"M1375_ES_IGO_11769_1_S59_L002_R1_001"] = as.matrix(counts[,"M1375_ES_IGO_11769_1_S59_L002_R1_001"] + counts[,"M1375_ES_IGO_11769_1_S59_L003_R1_001"])
counts[,"M1376_ES_IGO_11769_4_S67_L002_R1_001"] = as.matrix(counts[,"M1376_ES_IGO_11769_4_S67_L002_R1_001"] + counts[,"M1376_ES_IGO_11769_4_S67_L003_R1_001"])
counts[,"M1376_ES_IGO_11769_3_S66_L002_R1_001"] = as.matrix(counts[,"M1376_ES_IGO_11769_3_S66_L002_R1_001"] + counts[,"M1376_ES_IGO_11769_3_S66_L003_R1_001"])
counts[,"M1521_ES_IGO_11769_6_S69_L002_R1_001"] = as.matrix(counts[,"M1521_ES_IGO_11769_6_S69_L002_R1_001"] + counts[,"M1521_ES_IGO_11769_6_S69_L003_R1_001"])
counts[,"M1521_ES_IGO_11769_5_S68_L002_R1_001"] = as.matrix(counts[,"M1521_ES_IGO_11769_5_S68_L002_R1_001"] + counts[,"M1521_ES_IGO_11769_5_S68_L003_R1_001"])
counts[,"M1524_ES_IGO_11769_8_S71_L002_R1_001"] = as.matrix(counts[,"M1524_ES_IGO_11769_8_S71_L002_R1_001"] + counts[,"M1524_ES_IGO_11769_8_S71_L003_R1_001"])
counts[,"M1524_ES_IGO_11769_7_S70_L002_R1_001"] = as.matrix(counts[,"M1524_ES_IGO_11769_7_S70_L002_R1_001"] + counts[,"M1524_ES_IGO_11769_7_S70_L003_R1_001"])
counts[,"M2617_2mo_IGO_11769_12_S51_L002_R1_001"] = as.matrix(counts[,"M2617_2mo_IGO_11769_12_S51_L002_R1_001"] + counts[,"M2617_2mo_IGO_11769_12_S51_L003_R1_001"])
counts[,"M2617_2mo_IGO_11769_11_S50_L002_R1_001"] = as.matrix(counts[,"M2617_2mo_IGO_11769_11_S50_L002_R1_001"] + counts[,"M2617_2mo_IGO_11769_11_S50_L003_R1_001"])
counts[,"M2618_2mo_IGO_11769_10_S49_L002_R1_001"] = as.matrix(counts[,"M2618_2mo_IGO_11769_10_S49_L002_R1_001"] + counts[,"M2618_2mo_IGO_11769_10_S49_L003_R1_001"])
counts[,"M2618_2mo_IGO_11769_9_S72_L002_R1_001"] = as.matrix(counts[,"M2618_2mo_IGO_11769_9_S72_L002_R1_001"] + counts[,"M2618_2mo_IGO_11769_9_S72_L003_R1_001"])
counts[,"M2621_2mo_IGO_11769_20_S60_L002_R1_001"] = as.matrix(counts[,"M2621_2mo_IGO_11769_20_S60_L002_R1_001"] + counts[,"M2621_2mo_IGO_11769_20_S60_L003_R1_001"])
counts[,"M2621_2mo_IGO_11769_19_S58_L002_R1_001"] = as.matrix(counts[,"M2621_2mo_IGO_11769_19_S58_L002_R1_001"] + counts[,"M2621_2mo_IGO_11769_19_S58_L003_R1_001"])
counts[,"M2624_2mo_IGO_11769_16_S55_L002_R1_001"] = as.matrix(counts[,"M2624_2mo_IGO_11769_16_S55_L002_R1_001"] + counts[,"M2624_2mo_IGO_11769_16_S55_L003_R1_001"])
counts[,"M2624_2mo_IGO_11769_15_S54_L002_R1_001"] = as.matrix(counts[,"M2624_2mo_IGO_11769_15_S54_L002_R1_001"] + counts[,"M2624_2mo_IGO_11769_15_S54_L003_R1_001"])
counts[,"M2629_2mo_IGO_11769_14_S53_L002_R1_001"] = as.matrix(counts[,"M2629_2mo_IGO_11769_14_S53_L002_R1_001"] + counts[,"M2629_2mo_IGO_11769_14_S53_L003_R1_001"])
counts[,"M2629_2mo_IGO_11769_13_S52_L002_R1_001"] = as.matrix(counts[,"M2629_2mo_IGO_11769_13_S52_L002_R1_001"] + counts[,"M2629_2mo_IGO_11769_13_S52_L003_R1_001"])
counts[,"M2633_2mo_IGO_11769_18_S57_L002_R1_001"] = as.matrix(counts[,"M2633_2mo_IGO_11769_18_S57_L002_R1_001"] + counts[,"M2633_2mo_IGO_11769_18_S57_L003_R1_001"])
counts[,"M2633_2mo_IGO_11769_17_S56_L002_R1_001"] = as.matrix(counts[,"M2633_2mo_IGO_11769_17_S56_L002_R1_001"] + counts[,"M2633_2mo_IGO_11769_17_S56_L003_R1_001"])


# remove the technical replicates from the count matrix
counts <- counts[,-which(colnames(counts) %in% c('M1375_ES_IGO_11769_2_S65_L003_R1_001',
                                                 'M1375_ES_IGO_11769_1_S59_L003_R1_001',
                                                 'M1376_ES_IGO_11769_4_S67_L003_R1_001',
                                                 'M1376_ES_IGO_11769_3_S66_L003_R1_001',
                                                 'M1521_ES_IGO_11769_6_S69_L003_R1_001',
                                                 'M1521_ES_IGO_11769_5_S68_L003_R1_001',
                                                 'M1524_ES_IGO_11769_8_S71_L003_R1_001',
                                                 'M1524_ES_IGO_11769_7_S70_L003_R1_001',
                                                 'M2617_2mo_IGO_11769_12_S51_L003_R1_001',
                                                 'M2617_2mo_IGO_11769_11_S50_L003_R1_001',
                                                 'M2618_2mo_IGO_11769_10_S49_L003_R1_001',
                                                 'M2618_2mo_IGO_11769_9_S72_L003_R1_001',
                                                 'M2621_2mo_IGO_11769_20_S60_L003_R1_001',
                                                 'M2621_2mo_IGO_11769_19_S58_L003_R1_001',
                                                 'M2624_2mo_IGO_11769_16_S55_L003_R1_001',
                                                 'M2624_2mo_IGO_11769_15_S54_L003_R1_001',
                                                 'M2624_2mo_IGO_11769_15_S54_L003_R1_001',
                                                 'M2629_2mo_IGO_11769_14_S53_L003_R1_001',
                                                 'M2629_2mo_IGO_11769_13_S52_L003_R1_001',
                                                 'M2633_2mo_IGO_11769_18_S57_L003_R1_001',
                                                 'M2633_2mo_IGO_11769_17_S56_L003_R1_001'))]

colnames(counts)
dim(counts)

# remove the technical replicates from the phenotable
pheno <- pheno[-which(pheno$sample_names %in% c('M1375_ES_IGO_11769_2_S65_L003_R1_001',
                                                'M1375_ES_IGO_11769_1_S59_L003_R1_001',
                                                'M1376_ES_IGO_11769_4_S67_L003_R1_001',
                                                'M1376_ES_IGO_11769_3_S66_L003_R1_001',
                                                'M1521_ES_IGO_11769_6_S69_L003_R1_001',
                                                'M1521_ES_IGO_11769_5_S68_L003_R1_001',
                                                'M1524_ES_IGO_11769_8_S71_L003_R1_001',
                                                'M1524_ES_IGO_11769_7_S70_L003_R1_001',
                                                'M2617_2mo_IGO_11769_12_S51_L003_R1_001',
                                                'M2617_2mo_IGO_11769_11_S50_L003_R1_001',
                                                'M2618_2mo_IGO_11769_10_S49_L003_R1_001',
                                                'M2618_2mo_IGO_11769_9_S72_L003_R1_001',
                                                'M2621_2mo_IGO_11769_20_S60_L003_R1_001',
                                                'M2621_2mo_IGO_11769_19_S58_L003_R1_001',
                                                'M2624_2mo_IGO_11769_16_S55_L003_R1_001',
                                                'M2624_2mo_IGO_11769_15_S54_L003_R1_001',
                                                'M2624_2mo_IGO_11769_15_S54_L003_R1_001',
                                                'M2629_2mo_IGO_11769_14_S53_L003_R1_001',
                                                'M2629_2mo_IGO_11769_13_S52_L003_R1_001',
                                                'M2633_2mo_IGO_11769_18_S57_L003_R1_001',
                                                'M2633_2mo_IGO_11769_17_S56_L003_R1_001')),]


##---------------------------------------------
# Data importing and normalization
##---------------------------------------------
#make a DGE list
d <- DGEList(counts = counts, samples = pheno, group = pheno$condition)
d$samples$group <- relevel(x = d$samples$group, ref = "Ctrl")
d$samples$B_region <- as.factor(d$samples$B_region)
d$samples$conRegion <- as.factor(paste(d$samples$group, d$samples$B_region, sep = '.'))
d$samples$conRegionAge <- as.factor(paste(d$samples$group, d$samples$B_region, d$samples$age, sep = '.'))
d$samples$conAge <- as.factor(paste(d$samples$group, d$samples$age, sep = '.'))


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
pheno_2mo <- pheno[which(pheno$age %in% 'ES'),]
new_lcpm <- lcpm[,which(colnames(lcpm) %in% pheno_2mo$sample_names)]

df_pca <- prcomp(t(new_lcpm))
print(summary(df_pca))
df_out <- as.data.frame(df_pca$x)
df_out$sampleConditions <- pheno_2mo$condition
df_out$age <- pheno_2mo$age
df_out$B_region <- pheno_2mo$B_region


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
p <- p + geom_point(size = 6, stroke = 0.8)+ theme + scale_color_manual(values = c('#1339F4', '#C70039')) #+ xlab(percentage[1]) + ylab(percentage[2]) 
p 


##--------------------------------------
## Differential expression analysis 
##--------------------------------------
#combination__age <- factor(d$samples$conAge)
#design <- model.matrix(~0 + combination__age)
#colnames(design) <- levels(combination__age)


combination_B_region_age <- factor(d$samples$conRegionAge)
design <- model.matrix(~0 + combination_B_region_age)
colnames(design) <- levels(combination_B_region_age)

#group <- factor(d$samples$group)
#design <- model.matrix(~0 + group)
#colnames(design) <- levels(group)

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
    BS_2mo = VE.BS.2mo-Ctrl.BS.2mo,
    BS_ES = VE.BS.ES-Ctrl.BS.ES,
    Cort_2mo = VE.Cort.2mo-Ctrl.Cort.2mo,
    Cort_ES = VE.Cort.ES-Ctrl.Cort.ES,
    #CvsVE = VE-Ctrl,
    #CvsVE_BS = VE.BS-Ctrl.BS,
    #CvsVE_Cort = VE.Cort-Ctrl.Cort,
    levels = design)


fit <- glmQLFit(d, design, robust = T)
qlf <- glmQLFTest(fit, contrast=cont.matrix[,"Cort_ES"])
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

allGenes <- mapIds(org.Mm.eg.db, keys = rownames(sigGenes), keytype = "ENSEMBL", column="SYMBOL")
allGenes <- as.data.frame(allGenes)
testVlookup <- merge(sigGenes, allGenes, by.x="row.names", by.y="row.names") 

#write.table(x = testVlookup, file = "~/Desktop/UP_genes_BS_ES.txt", quote = F, sep = "\t", row.names = T)
#write.table(x = testVlookup, file = "~/Desktop/DownGenes_BS_ES.txt", quote = F, sep = "\t", row.names = T)

##-----------------------------------
## Progeny 
##-----------------------------------
#Gene names
gene_names <- mapIds(org.Mm.eg.db, keys = rownames(lcpm), keytype = "ENSEMBL", column="SYMBOL")
gene_names <- as.data.frame(gene_names)
#testVlookup <- merge(matrix, gene_names, by.x="row.names", by.y="row.names") 
#rownames(lcpm) <- gene_names$gene_names

matched = match(rownames(lcpm) , rownames(gene_names))
rownames(lcpm) = gene_names$gene_names[matched]

phenotable <- pheno %>% arrange(condition, B_region) %>% filter(age == "2mo", B_region == 'BS')
matrix_ordered <- lcpm[,phenotable$sample_names]


library(progeny)
pathways = progeny(matrix_ordered, scale=FALSE, organism = 'Mouse')
controls = phenotable$condition == "Ctrl"
ctl_mean = apply(pathways[controls,], 2, mean)
ctl_sd = apply(pathways[controls,], 2, sd)
pathways = t(apply(pathways, 1, function(x) x - ctl_mean))
pathways = apply(pathways, 1, function(x) x / ctl_sd)

library(dplyr)
result = apply(pathways, 1, function(x) {
    broom::tidy(lm(x ~ !controls)) %>%
        filter(term == "!controlsTRUE") %>%
        dplyr::select(-term)
})
mutate(bind_rows(result), pathway=names(result))

library(pheatmap)
myColor = colorRampPalette(c("Darkblue", "white","red"))(100)
pheatmap(pathways,fontsize=14, show_rownames = TRUE,
         color=myColor, main = "PROGENy", angle_col = 45, treeheight_col = 0,  
         border_color = NA)

df <- data.frame(Group = phenotable$condition)
rownames(df) <- colnames(matrix_ordered)
names(df) <- c( 'Group') 

Group        <-  c('#1339F4', '#C70039')
names(Group) <- c( "Ctrl", "VE")
anno_colors <- list(Group = Group)

pathways_selected <- pathways[c("EGFR", "JAK-STAT", "NFkB", "TNFa", "VEGF", "WNT"),]

pheatmap(pathways_selected, annotation_col = df, annotation_colors = anno_colors,  scale = "row", 
         show_rownames = T, 
         show_colnames = F, cluster_cols = F, cluster_rows = F)
##-----------------------------------
## Volcano plot 
##-----------------------------------
tt$table$threshold <- "NO"
tt$table$threshold[tt$table$logFC > 1.5 & tt$table$FDR < 0.05] <- "UP"
tt$table$threshold[tt$table$logFC <= -1.5 & tt$table$FDR < 0.05] <- "DOWN"
#tt$table$threshold = as.factor(abs(tt$table$logFC) >= 1.5 & tt$table$PValue < 3.721154e-06)

## Sort by ordered padj
allGenes <- mapIds(org.Mm.eg.db, keys = rownames(tt$table), keytype = "ENSEMBL", column="SYMBOL")
allGenes <- as.data.frame(allGenes)
table_ordered <- merge(tt$table, allGenes, by.x="row.names", by.y="row.names") 

# prder and label the most up, down and significant genes
table_ordered <- table_ordered %>% arrange(desc(logFC))
## Create a column to indicate which genes to label
table_ordered$genelabels <- ""
table_ordered$genelabels[1:30] <- table_ordered$allGenes[1:30]
table_ordered$genelabels <- ifelse(table_ordered$threshold == 'NO', "",table_ordered$genelabels)

table_ordered <- table_ordered %>% arrange(logFC)
table_ordered$genelabels[1:30] <- table_ordered$allGenes[1:30]
table_ordered$genelabels <- ifelse(table_ordered$threshold == 'NO', "",table_ordered$genelabels)

table_ordered <- table_ordered %>% arrange(FDR)
table_ordered$genelabels[1:30] <- table_ordered$allGenes[1:30]
table_ordered$genelabels <- ifelse(table_ordered$threshold == 'NO', "",table_ordered$genelabels)



g = ggplot(data = as.data.frame(table_ordered), aes(x=logFC, y=-log10(FDR), colour=threshold)) +
    geom_point(alpha=0.5, size=2)+
    #geom_text_repel(aes(x = logFC, y = -log10(FDR), label = genelabels)) +
    theme(axis.title.x=element_text(size=25), legend.position = "none",
          axis.title.y=element_text(size=25)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(plot.title = element_text(face="bold", size=20)) +
    theme(axis.text=element_text(size=16)) +
    theme(text = element_text(size=16)) +
    theme(axis.line = element_line(size = 0.5, colour = "black"))+
    theme(panel.border = element_rect(colour = "white",size = 0.1)) + 
    scale_color_manual(values=c('blue',"grey", "red")) +
    xlab("log2 fold change")  + ylab("-log10 FDR")
g



##-----------------------------------
## Neuroinflammation analysis
##-----------------------------------

genes_down <- c('PPP1R17',
                'CALB1',
                'PCP2',
                'GNG13',
                'SLC1A6',
                'RGS8',
                'NPAS4',
                'PVALB',
                'FOLR1',
                'NEUROG2',
                'EGR4',
                'GRID2IP',
                'PCP4',
                'TMEM200B',
                'CACNA1G',
                'C1QTNF4',
                'SCN4A',
                'GPR83',
                'CA8',
                'GRID2',
                'DUSP4')

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

sig_up = testVlookup[which(toupper(testVlookup$allGenes) %in% genes_up),]
sig_up
testVlookup[which(toupper(testVlookup$allGenes) %in% genes_down),]


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

phenotable_ordered <- pheno %>% arrange(condition, B_region) %>% filter(age == "2mo")
matrix_ordered <- matrix[,phenotable_ordered$sample_names]

gene_order <- c('Cxcl10', 'Il1b', 'Msr1', 'Socs3', 'Ctss', 'Pla2g4a','Il1a', 'C1qb', 'Ncf1','Tyrobp', 'Cx3cr1', 'Trem2', 'Csf1r')
matrix_ordered <- matrix_ordered[gene_order,]

pheatmap(matrix_ordered, annotation_col = df, annotation_colors = anno_colors,  scale = "row", 
         show_rownames = T, 
         show_colnames = F, cluster_cols = F, cluster_rows = F)

##-----------------------------------
## Heatmap of significant genes 
##-----------------------------------
matrix <- lcpm[which(rownames(lcpm) %in% rownames(SigGenesUpDown)),] #sig_up$allGene
#matrix <- lcpm[which(toupper(rownames(lcpm)) %in% toupper(genes_up)),] #genes_up
dim(matrix)

df <- data.frame(B_region = d$samples$B_region, Group = d$samples$group)
rownames(df) <- colnames(lcpm)
names(df) <- c("B_region", 'Group') 

Group        <-  c('#1339F4', '#C70039')
names(Group) <- levels(df$Group)

B_region        <- c('#22BE6E', '#FFC300')
names(B_region) <- levels(df$B_region)
anno_colors <- list(B_region = B_region, Group = Group)

phenotable_ordered <- pheno %>% arrange(condition, B_region) %>% filter(age == "2mo")
matrix_ordered <- matrix[,phenotable_ordered$sample_names]

#Gene names
gene_names <- mapIds(org.Mm.eg.db, keys = rownames(matrix_ordered), keytype = "ENSEMBL", column="SYMBOL")
gene_names <- as.data.frame(gene_names)
#testVlookup <- merge(matrix, gene_names, by.x="row.names", by.y="row.names") 
rownames(matrix_ordered) <- gene_names$gene_names


pheatmap(matrix_ordered, annotation_col = df, annotation_colors = anno_colors,  scale = "row", 
         show_rownames = T, 
         show_colnames = F, cluster_cols = F)



##---------------------------------------------
## Keren-Shaul et al DAM signature enrichment 
##---------------------------------------------

dam_sig <- read_table(file = "input/full_dam_signature.txt")
gene_names <- mapIds(org.Mm.eg.db, keys = rownames(lcpm), keytype = "ENSEMBL", column="SYMBOL")
gene_names <- as.data.frame(gene_names)
#testVlookup <- merge(lcpm, gene_names, by.x="row.names", by.y="row.names") 
rownames(lcpm) <- gene_names$gene_names

matrix <- lcpm[which(rownames(lcpm) %in% dam_sig$Gene_name),] #dam genes
dim(matrix)

data_dam = data.frame(t(matrix), group = d$samples$group, B_region = d$samples$B_region, age = d$samples$age)
melted = melt(data = data_dam, id.vars = c("group", "B_region", 'age'))
data_BS_2mo <- melted %>% filter(B_region == "BS" & age == "2mo")
data_CX_2mo <- melted %>% filter(B_region == "Cort" & age == "2mo")

data_BS_ES <- melted %>% filter(B_region == "BS" & age == "ES")
data_CX_ES <- melted %>% filter(B_region == "Cort" & age == "ES")


theme <- theme(panel.background = element_blank(),panel.border = element_rect(fill = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               strip.background=element_blank(),axis.text.x=element_text(colour="black"),
               axis.text.y = element_text(colour="black"), axis.text = element_text(size = 16), axis.ticks = element_line(colour = "black"),
               plot.margin = unit(c(1,1,1,1),"line"), legend.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 16, face = "bold")
               ,axis.title.y = element_text(color = "black", size = 16, face = "bold")) 
theme_update(plot.title = element_text(hjust = 0.5, face = "bold"))


ggplot(data = data_CX_ES, aes(x = group, y = value, fill = group)) + geom_boxplot() + theme + ylab("DAM signature") + scale_fill_manual(values = c('#2983d5', '#cb5840'))


##-----------------------------------
## Heatmap of DAM genes 
##-----------------------------------
# have checked ortho, LYZ2 doenst have a human ortho
dam1_genes <- c("CX3CR1", "P2RY12", "TMEM119", "TYROBP", "CTSB", "CTSD", "APOE", "B2M", "FTH1", "LYZ2")
dam2_genes <- c( "TREM2", "AXL", "CST7",
                 "CTSL", "LPL", "CD9", "CSF1", "CCL6", "ITGAX", "CLEC7A", "LILRB4", "TIMP2")

# Check how many of the dam genes are significant
allGenes <- mapIds(org.Mm.eg.db, keys = rownames(table), keytype = "ENSEMBL", column="SYMBOL")
allGenes <- as.data.frame(allGenes)
testVlookup_sigGenes <- merge(table, allGenes, by.x="row.names", by.y="row.names") 

common <- testVlookup_sigGenes[which(toupper(testVlookup_sigGenes$allGenes) %in% dam2_genes),]
common[order(common$FDR),]

#Gene names
gene_names <- mapIds(org.Mm.eg.db, keys = rownames(lcpm), keytype = "ENSEMBL", column="SYMBOL")
gene_names <- as.data.frame(gene_names)
rownames(lcpm) <- gene_names$gene_names
matrix <- lcpm[which(toupper(rownames(lcpm)) %in% dam2_genes),] #sig_up$allGene
#matrix <- lcpm[which(toupper(rownames(lcpm)) %in% toupper(genes_up)),] #genes_up
dim(matrix)

df <- data.frame(Age = d$samples$age, B_region = d$samples$B_region, Group = d$samples$group)
rownames(df) <- colnames(lcpm)
names(df) <- c('Age', "B_region", 'Group') 

Group        <-  c('#1339F4', '#C70039')
names(Group) <- levels(df$Group)

B_region        <- c('#22BE6E', '#FFC300')
names(B_region) <- levels(df$B_region)

Age       <- c('#282409', '#8ED64B')
names(Age) <- levels(as.factor(df$Age))

anno_colors <- list( Age = Age, B_region = B_region, Group = Group)

phenotable_ordered <- pheno %>% arrange(condition, B_region) %>% filter(age == "2mo")
matrix_ordered <- matrix[,phenotable_ordered$sample_names]

gene_order <- common[order(common$FDR),]$allGenes
matrix_ordered <- matrix_ordered[match(gene_order, rownames(matrix_ordered)), ]        # Reorder data frame

pheatmap(matrix_ordered, annotation_col = df, annotation_colors = anno_colors,  scale = "row", 
         show_rownames = T, 
         show_colnames = F, cluster_cols = F, cluster_rows = F)

# ##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ## GSEA for DAM genes 
# ##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#rnk.file <- read.csv(file = "Results/GSEA/BS_2mo_ranked_list.rnk", header = T, sep = '\t')
rnk.file <- read.csv(file = "Results/GSEA/BS_ES_ranked_list.rnk", header = T, sep = '\t')
dam_list <- list("dam1_genes" = dam1_genes, "dam2_genes" = dam2_genes)

# make it into a vector
ranks <- deframe(rnk.file)

# run GSEA
fgseaRes_DAM <- fgsea(pathways = dam_list, stats = ranks, minSize = 5, maxSize = 500)

fgseaResTidy_DAM <- fgseaRes_DAM %>%
    as_tibble() %>%
    arrange(desc(NES)) 

fgseaResTidy_filtered_DAM <- fgseaResTidy_DAM %>%
    as_tibble() %>% filter(padj <= 0.25)
fgseaResTidy_filtered_DAM

# plots the pathways
ggplot(fgseaResTidy_filtered_DAM, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=NES<1)) +  scale_fill_manual(values = c("#C70039", "#0024FF")) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Senescence signatures") + 
    theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plotEnrichment(dam_list[["dam2_genes"]],
               ranks) + labs(title="dam2_genes")

# ##------------------------------------------------------------------------------------------------------------------------------------------------------
# ## Heatmap of Senescent genes
# ##------------------------------------------------------------------------------------------------------------------------------------------------------
Segura_UP <- read.delim("/juno/work/geissmann/data/R_projects/wholeBrainHuman/Results/Senescence/Segura et al /Segura_UP.txt")
Segura_DOWN <- read.delim("/juno/work/geissmann/data/R_projects/wholeBrainHuman/Results/Senescence/Segura et al /Segura_DOWN.txt")
casella_UP <- read.csv("/juno/work/geissmann/data/R_projects/wholeBrainHuman/Results/Senescence/Cassella et al /casella_up.txt", sep="")
casella_DOWN <- read.delim("/juno/work/geissmann/data/R_projects/wholeBrainHuman/Results/Senescence/Cassella et al /casella_down.txt")
Fridman_UP <- read.delim("/juno/work/geissmann/data/R_projects/wholeBrainHuman/Results/Senescence/Fridman et al/Fridman_up.txt")
hu_UP <- read.delim("/juno/work/geissmann/data/R_projects/wholeBrainHuman/Results/Senescence/Hu et al /hu_up.txt")
Merad_UP <- read.csv("/juno/work/geissmann/data/R_projects/wholeBrainHuman/Results/Senescence/Merad/Merad_up.txt", sep="")

# Check how many of the dam genes are significant
allGenes <- mapIds(org.Mm.eg.db, keys = rownames(table), keytype = "ENSEMBL", column="SYMBOL")
allGenes <- as.data.frame(allGenes)
testVlookup_sigGenes <- merge(table, allGenes, by.x="row.names", by.y="row.names") 

common <- testVlookup_sigGenes[which(toupper(testVlookup_sigGenes$allGenes) %in% casella_UP$Gene),]
common[order(common$FDR),]

#Gene names
gene_names <- mapIds(org.Mm.eg.db, keys = rownames(lcpm), keytype = "ENSEMBL", column="SYMBOL")
gene_names <- as.data.frame(gene_names)
rownames(lcpm) <- gene_names$gene_names
matrix <- lcpm[which(toupper(rownames(lcpm)) %in% casella_UP$Gene),] #sig_up$allGene
dim(matrix)

df <- data.frame(Age = d$samples$age, B_region = d$samples$B_region, Group = d$samples$group)
rownames(df) <- colnames(lcpm)
names(df) <- c('Age', "B_region", 'Group') 

Group        <-  c('#1339F4', '#C70039')
names(Group) <- levels(df$Group)

B_region        <- c('#22BE6E', '#FFC300')
names(B_region) <- levels(df$B_region)

Age       <- c('#282409', '#8ED64B')
names(Age) <- levels(as.factor(df$Age))

anno_colors <- list( Age = Age, B_region = B_region, Group = Group)

phenotable_ordered <- pheno %>% arrange(condition, B_region) %>% filter(age == "ES")
matrix_ordered <- matrix[,phenotable_ordered$sample_names]

gene_order <- common[order(common$FDR),]$allGenes
matrix_ordered <- matrix_ordered[match(gene_order, rownames(matrix_ordered)), ]        # Reorder data frame

pheatmap(matrix_ordered, annotation_col = df, annotation_colors = anno_colors,  scale = "row", 
         show_rownames = T, 
         show_colnames = F, cluster_cols = F, cluster_rows = F)

# ##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ## GSEA for Senescence genes 
# ##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rnk.file <- read.csv(file = "Results/GSEA/BS_2mo_ranked_list.rnk", header = T, sep = '\t')
#rnk.file <- read.csv(file = "Results/GSEA/BS_ES_ranked_list.rnk", header = T, sep = '\t')
sen_list <- list("Segura_UP" = Segura_UP$Gene, "Segura_DOWN" = Segura_DOWN$Gene.name, "casella_UP" = casella_UP$Gene, "casella_DOWN" = casella_DOWN$Gene, 
                 "Fridman_UP" = Fridman_UP$Gene, "hu_UP" = hu_UP$Gene, "Merad_UP" = Merad_UP$Gene)

# make it into a vector
ranks <- deframe(rnk.file)

# run GSEA
fgseaRes_Sen <- fgsea(pathways = sen_list, stats = ranks, minSize = 5, maxSize = 500)

fgseaResTidy_Sen <- fgseaRes_Sen %>%
    as_tibble() %>%
    arrange(desc(NES)) 

fgseaResTidy_filtered_Sen <- fgseaResTidy_Sen %>%
    as_tibble() %>% filter(padj <= 0.25)
fgseaResTidy_filtered_Sen

# plots the pathways
ggplot(fgseaResTidy_filtered_Sen, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=NES<1)) +  scale_fill_manual(values = c("#C70039", "#0024FF")) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Senescence signatures") + 
    theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plotEnrichment(dam_list[["dam2_genes"]],
               ranks) + labs(title="dam2_genes")


##-------------------------------------------------
## Heatmap of common genes between Human and BS 2mo or BS ES 
##-------------------------------------------------
UP_human_mouse_common_2mo <- read.delim("/juno/work/geissmann/data/R_projects/bulk_microglia/Results/Human_Mouse_Common_Genes/UP_human_mouse_common_2mo.txt")
UP_human_mouse_common_ES <- read.delim("/juno/work/geissmann/data/R_projects/bulk_microglia/Results/Human_Mouse_Common_Genes/UP_human_mouse_common_ES.txt")

matrix <- lcpm[which(rownames(lcpm) %in% UP_human_mouse_common_ES$Mouse.Gene_ID),] #sig_up$allGene
dim(matrix)

df <- data.frame(Age = d$samples$age, B_region = d$samples$B_region, Group = d$samples$group)
rownames(df) <- colnames(lcpm)
names(df) <- c('Age', "B_region", 'Group') 

Group        <-  c('#1339F4', '#C70039')
names(Group) <- levels(df$Group)

B_region        <- c('#22BE6E', '#FFC300')
names(B_region) <- levels(df$B_region)

Age       <- c('#282409', '#8ED64B')
names(Age) <- levels(as.factor(df$Age))

anno_colors <- list( Age = Age, B_region = B_region, Group = Group)

phenotable_ordered <- pheno %>%  filter(age == "ES")  %>% arrange(condition, B_region) 
matrix_ordered <- matrix[,phenotable_ordered$sample_names]

#Gene names
gene_names <- mapIds(org.Mm.eg.db, keys = rownames(matrix_ordered), keytype = "ENSEMBL", column="SYMBOL")
gene_names <- as.data.frame(gene_names)
#testVlookup <- merge(matrix, gene_names, by.x="row.names", by.y="row.names") 
rownames(matrix_ordered) <- gene_names$gene_names

pheatmap(matrix_ordered, annotation_col = df, annotation_colors = anno_colors,  scale = "row", 
         show_rownames = T, 
         show_colnames = F, cluster_cols = F, cluster_rows = F)


##---------------------------------
## GSEA analysis pre-processing
##----------------------------------
library(biomaRt)
source(file = "/juno/work/geissmann/data/R_projects/bulk_microglia/source/create_GSEA_ranked_list.R")
source(file = "/juno/work/geissmann/data/R_projects/bulk_microglia/source/Mouse2Human.R")
# create ranked list for GSEA
allGenes <- mapIds(org.Mm.eg.db, keys = rownames(table), keytype = "ENSEMBL", column="SYMBOL")
allGenes <- as.data.frame(allGenes)
testVlookup <- merge(table, allGenes, by.x="row.names", by.y="row.names") 
testVlookup<- testVlookup %>% mutate_if(is.factor, as.character) %>% dplyr::filter(!is.na(allGenes))
ranks_RNAseq <- as.data.frame(create_GSEA_ranked_list(table = testVlookup))

## create the conversion table from mouse to human 
Mouse2HumanTable <- Mouse2Human(MouseGenes = ranks_RNAseq$GeneName)
Mouse2HumanTable <- Mouse2HumanTable[-which(Mouse2HumanTable$HGNC == ""),]
Mouse2HumanTable <- Mouse2HumanTable %>% distinct(HGNC, .keep_all = T)
mergeCols <- c("GeneName" = "MGI")
inner <- inner_join(ranks_RNAseq, Mouse2HumanTable, by = mergeCols)
ranked_list_human_orth <-  inner %>% dplyr::select(HGNC, rank)
#write.table(ranked_list_human_orth, file = "Results/DEA_results/BS_2mo_ranked_list.rnk", col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)
#write.table(ranked_list_human_orth, file = "Results/DEA_results/BS_ES_ranked_list.rnk", col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)

##---------------------------------
## fGSEA analysis
##----------------------------------

############################
## hallmarks pathway analysis 
############################
# load files 
rnk.file <- read.csv(file = "~/Desktop/BS_2mo_ranked_list.rnk", header = T, sep = '\t')
rnk.file_ES <- read.csv(file = "Results/DEA_results/BS_ES_ranked_list.rnk", header = T, sep = '\t') 
#rnk.file <- rnk.file_ES
pathways.hallmark <- gmtPathways("Results/GSEA/h.all.v7.4.symbols.gmt")

# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark %>% 
    head() %>% 
    lapply(head)

# make it into a vector
ranks <- deframe(rnk.file)

# run GSEA
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=15, maxSize=500)

fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) 

fgseaResTidy_filtered <- fgseaResTidy %>%
    as_tibble() %>% filter(padj <= 0.25) 
#write.table(fgseaResTidy_filtered[,1:7], file = "Results/GSEA/Hallmark_pathways.txt", col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)
#write.table(fgseaResTidy_filtered[,1:7], file = "Results/GSEA/BS_ESHallmark_pathways.txt", col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)


#fgseaResTidy %>% 
#    dplyr::select(-leadingEdge, -ES) %>% 
#    arrange(padj) %>%  DT::datatable()

# plots the pathways
ggplot(fgseaResTidy_filtered, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=NES<1)) + scale_fill_manual(values = c("#C70039", "#0024FF")) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Hallmark pathways") + 
    theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#hallmark_selected <- fgseaResTidy_filtered[c(1,2,5,7,11,12,15,18,20,21,22,24,27,30),]

#ggthemr('fresh')
# Selected plots the pathways
#ggplot(hallmark_selected, aes(reorder(pathway, NES), NES)) +
#    geom_col(aes(fill=padj<0.25)) +
#    coord_flip() +
#    labs(x="Pathway", y="Normalized Enrichment Score") + 
#    theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


plotEnrichment(pathways.hallmark[["HALLMARK_INFLAMMATORY_RESPONSE"]],
               ranks) + labs(title="HALLMARK_INFLAMMATORY_RESPONSE")


topPathwaysUp <- fgseaRes[NES > 0][head(order(padj), n=10), pathway]
topPathwaysDown <- fgseaRes[NES < 0][head(order(padj), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways.hallmark[topPathways], ranks, fgseaRes, gseaParam = 0.5)

# find the genes in the pathway that contribute to the enrichment score
fgseaResTidy[fgseaResTidy$pathway == "HALLMARK_INFLAMMATORY_RESPONSE", "leadingEdge"][[1]]

fgseaResTidy[fgseaResTidy$pathway == "HALLMARK_IL6_JAK_STAT3_SIGNALING", "leadingEdge"][[1]]


############################
## KEGG pathway analysis 
############################
pathways.KEGG <- gmtPathways("Results/GSEA/c2.cp.kegg.v7.4.symbols.gmt")
fgseaRes_KEGG <- fgsea(pathways=pathways.KEGG, stats=ranks, minSize=15, maxSize=500)

fgseaResTidy_filtered_KEGG <- fgseaRes_KEGG %>%
    as_tibble() %>% filter(padj <= 0.25) %>% arrange(padj)
#write.table(fgseaResTidy_filtered_KEGG[,1:7], file = "Results/GSEA/KEGG_pathways.txt", col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)
#write.table(fgseaResTidy_filtered_KEGG[,1:7], file = "Results/GSEA/BS_ES_KEGG_pathways.txt", col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)


ggplot(fgseaResTidy_filtered_KEGG, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=NES<1)) + scale_fill_manual(values = c("#C70039", "#0024FF")) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score", title = "KEGG pathways") + 
    theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


plotEnrichment(pathways.KEGG[["KEGG_ALZHEIMERS_DISEASE"]],ranks) + labs(title="KEGG_ALZHEIMERS_DISEASE")

# find the genes in the pathway that contribute to the enrichment score
fgseaResTidy_filtered_KEGG[fgseaResTidy_filtered_KEGG$pathway == "KEGG_ALZHEIMERS_DISEASE", "leadingEdge"][[1]]

############################
## GO pathway analysis 
############################

pathways.GO <- gmtPathways("Results/GSEA/c5.go.v7.4.symbols.gmt")
fgseaRes_GO <- fgsea(pathways=pathways.GO, stats=ranks, minSize=15, maxSize=500)

fgseaResTidy_filtered_GO <- fgseaRes_GO %>%
    as_tibble() %>% filter(padj <= 0.25) %>% arrange(padj)
#write.table(fgseaResTidy_filtered_GO[,1:7], file = "Results/GSEA/GO_pathways.txt", col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)
#write.table(fgseaResTidy_filtered_GO[,1:7], file = "Results/GSEA/BS_ES_GO_pathways.txt", col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)

ggplot(fgseaResTidy_filtered_GO[c(1,2,3,6,7,8,11:50),], aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=NES<1)) + scale_fill_manual(values = c("#C70039", "#0024FF")) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score", title = "Gene Ontology") + 
    theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


plotEnrichment(pathways.GO[["GOBP_CELLULAR_SENESCENCE"]], ranks) + labs(title="GOBP_CELLULAR_SENESCENCE")

# find the genes in the pathway that contribute to the enrichment score
fgseaResTidy_filtered_GO[fgseaResTidy_filtered_GO$pathway == "GOBP_CELLULAR_SENESCENCE", "leadingEdge"][[1]]


############################
## REACTOME pathway analysis 
############################

pathways.reac <- gmtPathways("Results/GSEA/c2.cp.reactome.v7.4.symbols.gmt")
fgseaRes_reac <- fgsea(pathways=pathways.reac, stats=ranks, minSize=15, maxSize=500)

fgseaResTidy_filtered_reac <- fgseaRes_reac %>%
    as_tibble() %>% filter(padj <= 0.25) %>% arrange(padj)
fgseaResTidy_filtered_reac
#write.table(fgseaResTidy_filtered_reac[,1:7], file = "Results/GSEA/REACTOME_pathways.txt", col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)
#write.table(fgseaResTidy_filtered_reac[,1:7], file = "Results/GSEA/BS_ES_REACTOME_pathways.txt", col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)

REAC_select <- fgseaResTidy_filtered_reac[c(5,9,10,12,21,32,55),]
REAC_select <- fgseaResTidy_filtered_reac[c(2,12,13,15,17,19,22,23,25,34,35,43,44,46,75,77),] 

#topPathwaysUp <- fgseaRes_reac[NES > 0][head(order(padj), n=30), c("pathway","pval","padj", "NES")]
#topPathwaysDown <- fgseaRes_reac[NES < 0][head(order(padj), n=5), c("pathway","pval","padj", "NES")]
#topPathways <- rbind(topPathwaysUp, topPathwaysDown)
#topPathways <- topPathways %>% as_tibble()

ggplot(REAC_select, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=NES<1)) + scale_fill_manual(values = c("#C70039", "#0024FF")) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score", title = "REACTOME") + 
    theme_light() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


plotEnrichment(pathways.reac[["REACTOME_CELLULAR_SENESCENCE"]], ranks) + labs(title="REACTOME_CELLULAR_SENESCENCE")

# find the genes in the pathway that contribute to the enrichment score
fgseaResTidy_filtered_reac[fgseaResTidy_filtered_reac$pathway == "REACTOME_CELLULAR_SENESCENCE", "leadingEdge"][[1]]

#save(fgseaRes, fgseaRes_KEGG, fgseaRes_GO, fgseaRes_reac, file = "Results/GSEA_results.RData")
save(fgseaRes, fgseaRes_KEGG, fgseaRes_GO, fgseaRes_reac, file = "Results/GSEA_results_for_BS_ES.RData")
#load("Results/GSEA/GSEA_results.RData")



# 
# ##-----------------------------------
# ## Gene overlap between conditions
# ##-----------------------------------
# library(GOplot)
# sigGenes_BS_ES <- testVlookup
# colnames(sigGenes_BS_ES)[7] <-c('ID')
# sigGenes_BS_2mo <- testVlookup
# colnames(sigGenes_BS_2mo)[7] <-c('ID')
# res_venn = GOVenn(sigGenes_BS_ES,sigGenes_BS_2mo, plot = F, label = c('BS_ES','BS_2mo'))
# venn_table= res_venn$table$AB
# venn_table_up = venn_table[venn_table$Trend == 'UP',]
# venn_table_down = venn_table[venn_table$Trend == 'DOWN',]
# 
# table_A_only <- res_venn$table$A_only
# table_B_only <- res_venn$table$B_only
# write.table(x = table_A_only, file = "~/Desktop/table_A_only_BS_ES.txt", quote = F, sep = "\t", row.names = T)
# write.table(x = table_B_only, file = "~/Desktop/table_B_only_BS_2mo.txt", quote = F, sep = "\t", row.names = T)
# 
# 
# ## Heatmap of commonly significant genes BS_ES and BS_2mo
# matrix <- lcpm[which(rownames(lcpm) %in% rownames(venn_table_down)),] #sig_up$allGene
# #matrix <- lcpm[which(toupper(rownames(lcpm)) %in% toupper(genes_up)),] #genes_up
# dim(matrix)
# 
# df <- data.frame(Age = d$samples$age, B_region = d$samples$B_region, Group = d$samples$group)
# rownames(df) <- colnames(lcpm)
# names(df) <- c('Age', "B_region", 'Group') 
# 
# Group        <-  c('#1339F4', '#C70039')
# names(Group) <- levels(df$Group)
# 
# B_region        <- c('#22BE6E', '#FFC300')
# names(B_region) <- levels(df$B_region)
# 
# Age       <- c('#282409', '#8ED64B')
# names(Age) <- levels(as.factor(df$Age))
# 
# anno_colors <- list( Age = Age, B_region = B_region, Group = Group)
# 
# 
# phenotable_ordered <- pheno %>% arrange(condition, B_region)
# matrix_ordered <- matrix[,phenotable_ordered$sample_names]
# 
# #Gene names
# gene_names <- mapIds(org.Mm.eg.db, keys = rownames(matrix_ordered), keytype = "ENSEMBL", column="SYMBOL")
# gene_names <- as.data.frame(gene_names)
# #testVlookup <- merge(matrix, gene_names, by.x="row.names", by.y="row.names") 
# rownames(matrix_ordered) <- gene_names$gene_names
# 
# 
# pheatmap(matrix_ordered, annotation_col = df, annotation_colors = anno_colors,  scale = "row", 
#          show_rownames = T, 
#          show_colnames = F, cluster_cols = F, cluster_rows = F)
# 

