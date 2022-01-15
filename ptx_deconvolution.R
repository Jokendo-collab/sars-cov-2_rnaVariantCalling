setwd("C:\\Users\\Javan\\Desktop\\sars-cov-2Project\\data\\human_sarscov_db")

suppressPackageStartupMessages(library(remotes))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplotify))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(hrbrthemes))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(ReactomePA))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(genefilter))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gage))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
#=========Celltypes deconvolution====================
# 1. Load the data
df = read.table("proteomeDeconvolution/CIBERSORTx_Job13_Adjusted.txt",header = T,sep = '\t')

# 2 Drop unwanted columns
df = select(df, -contains(c("RMSE","Correlation","P.value","Absolute.score..sig.score.")))

# 3 Visualize the data with the ggplot
ggplot(df, aes(x = Collection_site, y = Eosinophils, fill = Collection_site)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~Collection_site) +
  theme() +
  ggtitle("")

# 4. Plot the heatmap
dframe = read.table("proteomeDeconvolution/heatmap/CIBERSORTx_Job13_Adjusted.txt",header = T, row.names = 1,sep = '\t')

#Prepare metadata file
metadata = data.frame(
  row.names = colnames(dframe),
  Sample = c("BAL_1","BAL_2","BAL_3","BAL_4","BAL_5","Gargle_1","Gargle_2","Gargle_3","Naso_1",  
             "Naso_10","Naso_11","Naso_12","Naso_13","Naso_14","Naso_15","Naso_16","Naso_17","Naso_18" 
             ,"Naso_19","Naso_2","Naso_3","Naso_4","Naso_5","Naso_6","Naso_8","Naso_9","Urine_1", 
             "Urine_2","Urine_20","Urine_21","Urine_22","Urine_23","Urine_24","Urine_25","Urine_3","Urine_4" 
             ,"Urine_5","Urine_6","Urine_7","Urine_8","Urine_9"),
  Group = c("BAL","BAL","BAL","BAL","BAL","Gargle","Gargle","Gargle","Naso",  
            "Naso","Naso","Naso","Naso","Naso","Naso","Naso","Naso","Naso" 
            ,"Naso","Naso","Naso","Naso","Naso","Naso","Naso","Naso","Urine", 
            "Urine","Urine","Urine","Urine","Urine","Urine","Urine","Urine","Urine" 
            ,"Urine","Urine","Urine","Urine","Urine"))
#Convert the values to integers
dframe$BAL_1 = as.integer(dframe$BAL_1)
dframe$BAL_2 = as.integer(dframe$BAL_2)
dframe$BAL_3 = as.integer(dframe$BAL_3)
dframe$BAL_4 = as.integer(dframe$BAL_4)
dframe$BAL_5 = as.integer(dframe$BAL_5)
dframe$Gargle_1 = as.integer(dframe$Gargle_1)
dframe$Gargle_2 = as.integer(dframe$Gargle_2)
dframe$Gargle_3 = as.integer(dframe$Gargle_3)
dframe$Naso_1 = as.integer(dframe$Naso_1)
dframe$Naso_10 = as.integer(dframe$Naso_10)
dframe$Naso_11 = as.integer(dframe$Naso_11)
dframe$Naso_12 = as.integer(dframe$Naso_12)
dframe$Naso_13 = as.integer(dframe$Naso_13)
dframe$Naso_14 = as.integer(dframe$Naso_14)
dframe$Naso_15 = as.integer(dframe$Naso_15)
dframe$Naso_16 = as.integer(dframe$Naso_16)
dframe$Naso_17 = as.integer(dframe$Naso_17)
dframe$Naso_18 = as.integer(dframe$Naso_18)
dframe$Naso_19 = as.integer(dframe$Naso_19)
dframe$Naso_2 = as.integer(dframe$Naso_2)
dframe$Naso_3 = as.integer(dframe$Naso_3)
dframe$Naso_4 = as.integer(dframe$Naso_4)
dframe$Naso_5 = as.integer(dframe$Naso_5)
dframe$Naso_6 = as.integer(dframe$Naso_6)
dframe$Naso_8 = as.integer(dframe$Naso_8)
dframe$Naso_9 = as.integer(dframe$Naso_9)
dframe$Urine_1 = as.integer(dframe$Urine_1)
dframe$Urine_2 = as.integer(dframe$Urine_2)
dframe$Urine_20 = as.integer(dframe$Urine_20)
dframe$Urine_21 = as.integer(dframe$Urine_21)
dframe$Urine_22 = as.integer(dframe$Urine_22)
dframe$Urine_23 = as.integer(dframe$Urine_23)
dframe$Urine_24 = as.integer(dframe$Urine_24)
dframe$Urine_25 = as.integer(dframe$Urine_25)
dframe$Urine_3 = as.integer(dframe$Urine_3)
dframe$Urine_4 = as.integer(dframe$Urine_4)
dframe$Urine_5 = as.integer(dframe$Urine_5)
dframe$Urine_6 = as.integer(dframe$Urine_6)
dframe$Urine_7 = as.integer(dframe$Urine_7)
dframe$Urine_8 = as.integer(dframe$Urine_8)
dframe$Urine_9 = as.integer(dframe$Urine_9)


#Change the Sample and Group categories as factors
metadata$Sample = factor(metadata$Sample)
metadata$Group = factor(metadata$Group)

#Check column ID == number of raws in metadata
ncol(dframe) == nrow(metadata)

rownames(metadata) <- colnames(dframe) #match the rownames in metadata with colnames in the count matrix data

rownames(metadata) == colnames(dframe)



#Create Deseq object
dds1 <- DESeqDataSetFromMatrix(countData = dframe,
                               colData = metadata,
                               design = ~Group)
dds1 #dim: 255 entried and 12 variables 

#================trial================
kdf = read.table("proteomeDeconvolution/heatmap/CIBERSORTx_Job13_Adjusted.txt",row.names = 1,header = T,sep = '\t')

g3<-melt(kdf)
plot1<-ggplot(g3, aes(variable, Mixture, fill= value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text = element_text(size = 7)) +
  ylab("Variant types") +
  xlab("Variant frequency") +
  ggtitle(" ") +
  scale_fill_distiller(palette = 'Spectral')

plot1

#===============ComplexHeatmaps=============
Heatmap(kdf)

Heatmap(kdf, name = "Abundance",column_km = 3,column_dend_side = "top",clustering_distance_rows = "pearson",row_km = 2,
        column_title = " ",ht_opt(RESET = TRUE))

#============PCA=============
library("factoextra")
library("FactoMineR")
pca.data <- PCA(kdf[,-1], scale.unit = TRUE, graph = FALSE)

fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 70))

fviz_pca_var(pca.data, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE) 

#=========================================================
kdf2 = read.table("proteomeDeconvolution/heatmap/Adjusted.txt",row.names = 1,header = T,sep = '\t')

wdbc.pr <- prcomp(kdf2[c(2:23)], center = TRUE, scale = TRUE)
summary(wdbc.pr)

screeplot(wdbc.pr, type = "l", npcs = 15, main = "Screeplot of the first 10 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)

cumpro <- cumsum(wdbc.pr$sdev^2 / sum(wdbc.pr$sdev^2))
plot(cumpro[0:15], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 6, col="blue", lty=5)
abline(h = 0.88759, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC6"),
       col=c("blue"), lty=5, cex=0.6)

plot(wdbc.pr$x[,1],wdbc.pr$x[,2], xlab="PC1 (33%)", ylab = "PC2 (28%)", main = "PC1 / PC2 - plot")

fviz_pca_ind(wdbc.pr, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = kdf2$Collection_site, 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Body site") +
  ggtitle(" ") +
  theme(plot.title = element_text(hjust = 0.8))










