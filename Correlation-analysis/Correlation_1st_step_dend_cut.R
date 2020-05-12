




#######################################################################
#T cell (smaller k) correlation network analysis - NC v VC and NU v VU
#######################################################################
# install.packages("Hmisc")
# install.packages("dendextend")
# install.packages("RColorBrewer")
# install.packages("ppcor")
# install.packages("qgraph")
#  install.packages("flashClust")
# install.packages("preprocessCore")
# install.packages("dplyr")
# install.packages("corrplot")
# install.packages("stats")
# install.packages("readxl")
# install.packages("openxlsx")
# install.packages("pheatmap")
# install.packages("dplyr")


library(dplyr)
library(Hmisc)
library(ppcor)
library(corrplot)
library(qgraph)
library(stats)
library(readxl)
library(openxlsx)
library(flashClust)
library(pheatmap)
library(dendextend)
library(RColorBrewer)
library(grid)
library(gtable)



#------------------- Read Data ---------------------------------------

dataD10T <- read.xlsx("HistoScore and Phenotype Data_T Cell Clusters(smaller K)_DR3.xlsx",
                      sheet = "Day 10")
dataD10T <- na.omit(dataD10T)
clusterD10T <- dplyr::select(dataD10T, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD10T)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD10T)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD10T)
prefix = "T D10 Cluster"
suffix = seq(1:n)
colnames(clusterD10T)[1:ncol(clusterD10T)] = paste(prefix, suffix, sep = " ")


dataD35T <- read.xlsx("HistoScore and Phenotype Data_T Cell Clusters(smaller K)_DR3.xlsx",
                      sheet = "Day 35")
dataD35T = na.omit(dataD35T)
clusterD35T <- dplyr::select(dataD35T, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD35T)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD35T)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD35T)
prefix = "T D35 Cluster"
suffix = seq(1:n)
colnames(clusterD35T)[1:ncol(clusterD35T)] = paste(prefix, suffix, sep = " ")


dataD44T <- read.xlsx("HistoScore and Phenotype Data_T Cell Clusters(smaller K)_DR3.xlsx",
                      sheet = "Day 44")
clusterD44T <- dplyr::select(dataD44T, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD44T)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD44T)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD44T)
prefix = "T D44 Cluster"
suffix = seq(1:n)
colnames(clusterD44T)[1:ncol(clusterD44T)] = paste(prefix, suffix, sep = " ")

dataD51T <- read.xlsx("HistoScore and Phenotype Data_T Cell Clusters(smaller K)_DR3.xlsx",
                      sheet = "Day 51")
clusterD51T <- dplyr::select(dataD51T, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD51T)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD51T)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD51T)
prefix = "T D51 Cluster"
suffix = seq(1:n)
colnames(clusterD51T)[1:ncol(clusterD51T)] = paste(prefix, suffix, sep = " ")


HistoData = dplyr::select(dataD51T, "Spleen.Burdon.(copies/mg)",
                          "Spleen/Body.wt..(ratio)",
                          "HeartScore",
                          "LiverScore",
                          "SpleenScore", 
                          "LungScore", "Day10.Antibody.OD",
                          "Day24.Antibody.OD",
                          "Day35.Antibody.OD")

colnames(HistoData) = c("Spleen Burden (copies/mg)", 
                        "Spleen/Body wt. (ratio)",
                        "Heart score",
                        "Liver score",
                        "Spleen score",
                        "Lung score",
                        "Day10 Antibody OD",
                        "Day24 Antibody OD",
                        "Day35 Antibody OD")

prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(HistoData)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(HistoData)[17:31] <- paste(prefix, suffix, sep = " ")




dataD10B <- read.xlsx("HistoScore and Phenotype Data_B Cell Clusters(smallerK)_DR3.xlsx",
                      sheet = "Day 10")
dataD10B <- na.omit(dataD10B)
clusterD10B <- dplyr::select(dataD10B, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD10B)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD10B)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD10B)
prefix = "B D10 Cluster"
suffix = seq(1:n)
colnames(clusterD10B)[1:ncol(clusterD10B)] = paste(prefix, suffix, sep = " ")


dataD35B <- read.xlsx("HistoScore and Phenotype Data_B Cell Clusters(smallerK)_DR3.xlsx",
                      sheet = "Day 35")
dataD35B = na.omit(dataD35B)
clusterD35B <- dplyr::select(dataD35B, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD35B)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD35B)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD35B)
prefix = "B D35 Cluster"
suffix = seq(1:n)
colnames(clusterD35B)[1:ncol(clusterD35B)] = paste(prefix, suffix, sep = " ")


dataD44B <- read.xlsx("HistoScore and Phenotype Data_B Cell Clusters(smallerK)_DR3.xlsx",
                      sheet = "Day 44")
clusterD44B <- dplyr::select(dataD44B, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD44B)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD44B)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD44B)
prefix = "B D44 Cluster"
suffix = seq(1:n)
colnames(clusterD44B)[1:ncol(clusterD44B)] = paste(prefix, suffix, sep = " ")

dataD51B <- read.xlsx("HistoScore and Phenotype Data_B Cell Clusters(smallerK)_DR3.xlsx",
                      sheet = "Day 51")
clusterD51B <- dplyr::select(dataD51B, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD51B)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD51B)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD51B)
prefix = "B D51 Cluster"
suffix = seq(1:n)
colnames(clusterD51B)[1:ncol(clusterD51B)] = paste(prefix, suffix, sep = " ")


dataD10I <- read.xlsx("HistoScore and Phenotype Data_Innate Cell Clusters(mixed K)_DR3.xlsx",
                      sheet = "Day 10")
dataD10I <- na.omit(dataD10I)
clusterD10I <- dplyr::select(dataD10I, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD10I)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD10I)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD10I)
prefix = "I D10 Cluster"
suffix = seq(1:n)
colnames(clusterD10I)[1:ncol(clusterD10I)] = paste(prefix, suffix, sep = " ")


dataD35I <- read.xlsx("HistoScore and Phenotype Data_Innate Cell Clusters(mixed K)_DR3.xlsx",
                      sheet = "Day 35")
dataD35I = na.omit(dataD35I)
clusterD35I <- dplyr::select(dataD35I, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD35I)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD35I)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD35I)
prefix = "I D35 Cluster"
suffix = seq(1:n)
colnames(clusterD35I)[1:ncol(clusterD35I)] = paste(prefix, suffix, sep = " ")


dataD44I <- read.xlsx("HistoScore and Phenotype Data_Innate Cell Clusters(mixed K)_DR3.xlsx",
                      sheet = "Day 44")
clusterD44I <- dplyr::select(dataD44I, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD44I)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD44I)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD44I)
prefix = "I D44 Cluster"
suffix = seq(1:n)
colnames(clusterD44I)[1:ncol(clusterD44I)] = paste(prefix, suffix, sep = " ")

dataD51I <- read.xlsx("HistoScore and Phenotype Data_Innate Cell Clusters(mixed K)_DR3.xlsx",
                      sheet = "Day 51")
clusterD51I <- dplyr::select(dataD51I, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD51I)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD51I)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD51I)
prefix = "I D51 Cluster"
suffix = seq(1:n)
colnames(clusterD51I)[1:ncol(clusterD51I)] = paste(prefix, suffix, sep = " ")


Category = dplyr::select(dataD51T, starts_with("Category"))
ClusterData = cbind(clusterD10T, 
                    clusterD35T, 
                    clusterD44T, 
                    clusterD51T,
                    clusterD10B,
                    clusterD35B,
                    clusterD44B,
                    clusterD51B,
                    clusterD10I,
                    clusterD35I,
                    clusterD44I,
                    clusterD51I,
                    Category)
HistoData = cbind(HistoData, Category)

ClusterdataNCvVC =  subset(ClusterData, Category == "Naïve Challenge" | Category == "Vaccinated Challenge")
HistoDataNCvVC = subset(HistoData, Category == "Naïve Challenge" | Category == "Vaccinated Challenge")
ClusterdataNUvVU = subset(ClusterData, Category == "Naïve Unchallenge" | Category == "Vaccinated Unchallenge")
HistoDataNUvVU = subset(HistoData, Category == "Naïve Unchallenge" | Category == "Vaccinated Unchallenge")


# ------------- Test for Normality -----------------------
# 
# shapiro = apply(ClusterData[,-which(names(ClusterData)=="Category")], 2, shapiro.test)
# shapiro
# p.value = unlist(lapply(shapiro, function(x) x$p.value))
# p.value = as.data.frame(p.value)
# write.xlsx(p.value, "D10-D52_Shapiro-Wilks_AllClusters.xlsx", row.names = TRUE)



########################################################
########################################################
#NCvVC Correlation Analysis
########################################################
########################################################



#-------- Use Spearmans for Physical Correlations -------

ClusterdataNCvVC = as.matrix(ClusterdataNCvVC[,-which(names(ClusterdataNCvVC)=="Category")])
HistoDataNCvVC = as.matrix(HistoDataNCvVC[,-which(names(HistoDataNCvVC)=="Category")])
TotalCorr_NCvVC = rcorr(ClusterdataNCvVC, HistoDataNCvVC, type="spearman")
Pvalue_NCvVC = TotalCorr_NCvVC$P
Coeff_NCvVC = TotalCorr_NCvVC$r

Coeff_NCvVC = as.data.frame(Coeff_NCvVC)
Pvalue_NCvVC = as.data.frame(Pvalue_NCvVC)
write.xlsx(Coeff_NCvVC, "D10-D51_Coefficients NCvVC Matrix.xlsx", row.names = TRUE)
write.xlsx(Pvalue_NCvVC, "D10-D51_Correlation NCvVC P-values.xlsx", row.names = TRUE)

#-----------Extract Cluster Assignments-----------------------
input_data_clustering = dplyr::select(Coeff_NCvVC, "Spleen Burden (copies/mg)",
                                      "Spleen/Body wt. (ratio)",
                                      "Heart score",
                                      "Liver score",
                                      "Spleen score",
                                      "Lung score",
                                      "Day10 Antibody OD",
                                      "Day24 Antibody OD",
                                      "Day35 Antibody OD")
number_clusters = 228

input_data_clustering_final = as.matrix(input_data_clustering[1:number_clusters,])

#Reorder names to match prism graphs
reorder.names = c("Day10 Antibody OD",
                  "Day24 Antibody OD",
                  "Day35 Antibody OD",
                  "Spleen Burden (copies/mg)",
                  "Spleen/Body wt. (ratio)",
                  "Spleen score",
                  "Lung score",
                  "Liver score",
                  "Heart score")
input_data_clustering_final = input_data_clustering_final[,reorder.names]

#Rename to match prism graphs
temp_colnames = c("Antibody D10",
                  "Antibody D24",
                  "Antibody D35",
                  "Spleen Burden",
                  "Spleen/Body Wgt",
                  "Spleen Score",
                  "Lung Score",
                  "Liver Score",
                  "Heart Score")
colnames(input_data_clustering_final) <- temp_colnames


#Plot heatmap
source("modified_pheatmap_extra_space.R")
# library(dplyr)
# input_data_clustering_final <- read.xlsx("correlation_matrix_data.xlsx", sheet = 1, rowNames = T)
#input_data_clustering_final <- as.matrix(input_data_clustering_final)

#input_data_clustering_final$Heart.Score
jpeg(paste("D10-D51_correlation_heatmap_7modules_ExtraSpace.jpeg"),units = "cm",res = 300,width = 50,height = 90)
pheatmap_modified(
  
  #dplyr::select(input_data_clustering_final, -c("Count", "IFNg", "CD184", "CD14", "IL17a")),
  input_data_clustering_final,
  
  
  scale = "none",
  
  color = colorRampPalette(c("#663366", "white", "#990000"))(100),
  
  border_color = "Black",
  
  
  
  clustering_distance_rows = "euclidean",
  
  # clustering_distance_cols = "euclidean",
  cluster_cols = F,
  
  clustering_method = "ward.D2",
  
  
  
  cellwidth = 25,
  
  cellheight = 10,
  
  
  
  treeheight_row = 200,
  
  treeheight_col = 50,
  
  
  
  fontsize = 9.5,
  
  
  
  cutree_rows = 2.7,
  
  display_numbers = F,
  number_color = "white"
  
  
  #width = 800
  
)

dev.off()





















