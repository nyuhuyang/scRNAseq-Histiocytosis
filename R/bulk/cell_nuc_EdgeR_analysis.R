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

pheno$protocol <- 'cell'


# ##---------------------------------------------
# # Load data and make Phenotable for Nuc
# ##---------------------------------------------

all_counts_nuc <- fc$counts
colnames(all_counts)

counts_nuc <- all_counts[,41:48]

list <- lapply(X = colnames(counts_nuc), FUN = strsplit, split = "[.]" )
list <- lapply(list, unlist)

# rename colnames
for(i in 1:length(list)){
    colnames(counts_nuc)[i] <- paste(list[[i]][11], list[[i]][12], list[[i]][13], sep = '_') #sometimes 10
}

# make pheno table
pheno_nuc <- data.frame(sample_names = colnames(counts_nuc))
for(i in 1:length(list)){
    pheno_nuc$condition[i] <- list[[i]][12]
}


for(i in 1:length(list)){
    print(list[[i]][13])
    list2 <- lapply(X = list[[i]][13], FUN = strsplit, split = "[_]" )
    list2 <- lapply(list2, unlist)
    pheno_nuc$B_region[i] <- list2[[1]][1]
}

#pheno$age <- ifelse(grepl(pattern = "2mo", x = pheno$sample_names), "2mo", "ES")
pheno_nuc$OG_name <- colnames(all_counts[,41:48])



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
pheno_nuc <- pheno_nuc[-which(pheno_nuc$sample_names %in% c('M2617_VE_BS_IGO_11769_24_S64_L003_R1_001',
                                                'M2617_VE_Cort_IGO_11769_23_S63_L003_R1_001',
                                                'M2618_Ctrl_BS_IGO_11769_22_S62_L003_R1_001',
                                                'M2618_Ctrl_Cort_IGO_11769_21_S61_L003_R1_001')),]
pheno_nuc$age <- '2mo'
pheno_nuc$protocol <- 'nuc'
pheno_nuc <- pheno_nuc[,c('sample_names', 'condition' ,'B_region', 'age', 'OG_name')] 

pheno_all <- rbind(pheno, pheno_nuc)
counts_all <- cbind(counts, counts_nuc)

##---------------------------------------------
# Data importing and normalization
##---------------------------------------------
#make a DGE list
d <- DGEList(counts = counts_all, samples = pheno_all, group = pheno_all$condition)
d$samples$group <- relevel(x = d$samples$group, ref = "Ctrl")
d$samples$B_region <- as.factor(d$samples$B_region)
d$samples$conRegion <- as.factor(paste(d$samples$group, d$samples$B_region, sep = '.'))
d$samples$conRegionAge <- as.factor(paste(d$samples$group, d$samples$B_region, d$samples$age, sep = '.'))
d$samples$conAge <- as.factor(paste(d$samples$group, d$samples$age, sep = '.'))
d$samples$protocol <- pheno_all$protocol

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
g <- ggplot(data = d$samples, aes(x = sample_names, y = lib.size, fill = protocol)) + geom_bar(stat="identity") 
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

pheno_2mo <- pheno_all[which(pheno_all$age %in% '2mo'),]
new_lcpm <- lcpm[,pheno_2mo$sample_names]

df_pca <- prcomp(t(new_lcpm))
print(summary(df_pca))
df_out <- as.data.frame(df_pca$x)
df_out$sampleConditions <- pheno_2mo$condition
df_out$age <- pheno_2mo$age
df_out$B_region <- pheno_2mo$B_region
df_out$protocol <- pheno_2mo$protocol


theme <- theme(panel.background = element_blank(),panel.border = element_rect(fill = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               strip.background=element_blank(),axis.text.x=element_text(colour="black"),
               axis.text.y = element_text(colour="black"), axis.text = element_text(size = 12), axis.ticks = element_line(colour = "black"),
               plot.margin = unit(c(1,1,1,1),"line"), legend.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 12, face = "bold")
               ,axis.title.y = element_text(color = "black", size = 12, face = "bold")) 
theme_update(plot.title = element_text(hjust = 0.5, face = "bold"))


p <- ggplot(df_out,aes(x = PC1, y = PC2, color = protocol))
p <- p + geom_point(size = 3, stroke = 0.8)+ theme + scale_color_manual(values = c('#00FF4A', '#FF00B5')) #+ xlab(percentage[1]) + ylab(percentage[2]) 
p 

