# install.packages("ggrepel")
# install.packages("Rtsne")

library(Rtsne)
library(openxlsx)
library(ggplot2)
library(ggrepel)
#library(ggsci)
##ploting function##
plot_tsne = function(data, type, clinical_information, label, label_tsne){
  
  # data = d_tsne_1
  # type = "type"
  # clinical_information = clinical
  # label = "cell.name"
  # label_tsne = label_tsne
  # 
  # colnames(d_tsne_1)
  
  plot = ggplot(data, aes_string(x="x", y="y", color = clinical_information, shape = type, size= "-log10(p.value)")) +
    geom_point() +
    scale_colour_gradient2(
      low="#0005a0",
      mid="#f7f7f7",
      high="#bf1400")+
    xlab("") + ylab("") +
    ggtitle("") +
    theme_light(base_size=15) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "vertical",
          legend.position = "right",
          legend.box = "vertical")
  if(label_tsne == TRUE){
    plot = plot +
      ggrepel::geom_text_repel(data = data,
                               ggplot2::aes_string(label = label),
                               size          = 5,
                               box.padding   = grid::unit(0.2, "lines"),
                               point.padding = grid::unit(0.2, "lines"))
  }
  return(plot)
}
##data_dealing function##
data_dealing <- function(Coefficients.matrix,
                         p.value.matrix, 
                         clinical,
                         grouping.file,
                         #p.value = 0.05,
                         label_tsne = TRUE,
                         metacluster = 5,
                         label_metacluster = 10){
  
# Coefficients.matrix = "D10-D51_Coefficients NCvVC Matrix.xlsx"
# p.value.matrix = "D10-D51_Correlation NCvVC P-values.xlsx"
# clinical = "Spleen.Burden.copies.mg"
# grouping.file = "group_data.xlsx"
# # p.value = 0.06,
# label_tsne = FALSE
# metacluster = 5
# label_metacluster = 3
  
  All.data <- read.xlsx(Coefficients.matrix,
                        sheet = 1,
                        rowNames = T)
  All.P.value <- read.xlsx(p.value.matrix,
                           sheet = 1,
                           rowNames = T)
  group.data <- read.xlsx(grouping.file,
                          sheet = 1,
                          rowNames = F)
  
  #rownames(All.data) == rownames(group.data)
  #All.data <- cbind(All.data, group.data)
  
  Cell.data <- All.data[which(substr(rownames(All.data), 1, 1) %in% c("T", "B", "I")), ]
  Cell.data.with.clinical <- Cell.data
  Cell.data <- Cell.data[, which(substr(colnames(Cell.data), 1, 1) %in% c("T", "B", "I"))]
  
  P.value.data <- All.P.value[which(substr(rownames(All.P.value), 1, 1) %in% c("T", "B", "I")), ]
  P.value.data.with.clinical <- P.value.data
  
  colnames(group.data)[1:2] <- c("cell.name", "group")
  group.data <- group.data[which(substr(group.data$cell.name, 1, 1) %in% c("T", "B", "I")), ]
  group.data = dplyr::select(group.data, c("cell.name","group"))
  group.data$group <- gsub("group_", "Module ", group.data$group, fixed = T)
  group.data$cell.name <- gsub(" ", ".", group.data$cell.name, fixed = T)
  #rownames(group.data) <- group.data$cluster
  #group.data <- group.data[, -1]
  
  for(i in 1: nrow(Cell.data.with.clinical)){
    #i =1
    if(substr(rownames(Cell.data.with.clinical)[i], 1, 1) == "T"){
      Cell.data.with.clinical$type[i] = "T cell"
    }else if(substr(rownames(Cell.data.with.clinical)[i], 1, 1) == "B"){
      Cell.data.with.clinical$type[i] = "B cell"
    }else if(substr(rownames(Cell.data.with.clinical)[i], 1, 1) == "I"){
      Cell.data.with.clinical$type[i] = "Innate cell"
    }
  }
  
  #colnames(Cell.data.with.clinical)
  #Cell.data.with.clinical$type
  
  
  ###########################################################
  ## Notice
  ## The dataset "Cell.data.with.clinical" should not be used directly 
  ## for tsen as it catained a column named type which is string variant.
  ## Logically, a string variant should not be used for data dimentional
  ## reduction, although this code doesn't produce any error.
  
  ###########################################################
  
  Cell.data.tsne <- dplyr::select(Cell.data.with.clinical, c("Spleen.Burden.copies.mg",
                                                             "Spleen.Body.wt.ratio",
                                                             "Heart.score",
                                                            "Liver.score",
                                                            "Spleen.score",
                                                            "Lung.score",
                                                            "Day10.Antibody.OD",
                                                            "Day24.Antibody.OD",
                                                            "Day35.Antibody.OD"))
  
  set.seed(1000)
  tsne <- Rtsne(#Cell.data.with.clinical,
    Cell.data.tsne, 
    dims = 2,
    perplexity=15,
    verbose=TRUE,
    max_iter = 10000,
    normalize = FALSE,
    pca=F, theta = 0, exaggeration_factor = 12)
  
  d_tsne_1 <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2])
  
  fit_cluster_kmeans=kmeans(scale(d_tsne_1), metacluster)
  fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))
  
  d_tsne_1$cl_kmeans = factor(fit_cluster_kmeans$cluster)
  d_tsne_1$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k = metacluster))
  d_tsne_1$type = factor(Cell.data.with.clinical$type, levels = unique(Cell.data.with.clinical$type))
  
  d_tsne_1[clinical] = Cell.data.with.clinical[clinical]
  d_tsne_1$cell.name = gsub(" ", ".", rownames(Cell.data.with.clinical), fixed = T)
  d_tsne_1$cell.name.full = d_tsne_1$cell.name
  
  d_tsne_1$p.value = P.value.data.with.clinical[[clinical]]
  d_tsne_1 = merge(d_tsne_1,group.data)
  
  for(i in 1: nrow(d_tsne_1)){
    #i =1
    if(d_tsne_1$cl_hierarchical[i] != label_metacluster){
      d_tsne_1$cell.name[i] = NA
    }
  }
  write.xlsx(d_tsne_1, paste(clinical, "_tsne_data.xlsx", sep = ""), row.names = F)
  #tiff("preliminery_view.tiff", res = 150)
  preliminery_view = ggplot(d_tsne_1, aes_string(x="x", y="y", color = "type", shape = "type", size= "-log10(p.value)")) +
    geom_point() +
    xlab("") + ylab("") +
    ggtitle("") +
    theme_light(base_size=15) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "vertical",
          legend.position = "right",
          legend.box = "vertical")
  ggsave("preliminery_view.jpeg", plot = preliminery_view, width = 10, height = 10)
  #dev.off()
  kmeans_view = ggplot(d_tsne_1, aes_string(x="x", y="y", color = "cl_kmeans", shape = "type", size= "-log10(p.value)")) +
    geom_point() +
    xlab("") + ylab("") +
    ggtitle("") +
    theme_light(base_size=15) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "vertical",
          legend.position = "right",
          legend.box = "vertical")
  ggsave("kmeans_view.jpeg", plot = kmeans_view, width = 10, height = 10)
  
  hierarchical_view = ggplot(d_tsne_1, aes_string(x="x", y="y", color = "cl_hierarchical", shape = "type", size= "-log10(p.value)")) +
    geom_point() +
    xlab("") + ylab("") +
    ggtitle("") +
    theme_light(base_size=15) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "vertical",
          legend.position = "right",
          legend.box = "vertical")
  ggsave("hierarchical_view.jpeg", plot = hierarchical_view, width = 10, height = 10)
  
  group_view = ggplot(d_tsne_1, aes_string(x="x", y="y", color = "group", shape = "type", size= 2)) +
    geom_point() +
    xlab("") + ylab("") +
    ggtitle("") +
    theme_light(base_size=15) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "vertical",
          legend.position = "right",
          legend.box = "vertical")
  ggsave("group_view.jpeg", plot = group_view, width = 10, height = 10)
  
  plot.clinical = plot_tsne(data = d_tsne_1,
                            type = "type",
                            clinical_information = clinical,
                            label = "cell.name",
                            label_tsne = label_tsne)
  plot.clinical
}
#?ggsave
data_dealing(Coefficients.matrix = "D10-D51_Coefficients NCvVC Matrix.xlsx",
             p.value.matrix = "D10-D51_Correlation NCvVC P-values.xlsx",
             clinical = "Spleen.Burden.copies.mg",
             grouping.file = "group_data.xlsx",
             # p.value = 0.06,
             label_tsne = FALSE,
             metacluster = 5,
             label_metacluster = 3)




#run loop to export all data tables
clinical_list = as.list(c("Spleen.Burden.copies.mg",
  "Spleen.Body.wt.ratio",
  "Heart.score",
  "Liver.score",
  "Spleen.score",
  "Lung.score",
  "Day10.Antibody.OD",
  "Day24.Antibody.OD",
  "Day35.Antibody.OD"))
for (i in clinical_list){
  data_dealing(Coefficients.matrix = "D10-D51_Coefficients NCvVC Matrix.xlsx",
               p.value.matrix = "D10-D51_Correlation NCvVC P-values.xlsx",
               clinical = i,
               grouping.file = "group_data.xlsx",
               # p.value = 0.06,
               label_tsne = FALSE,
               metacluster = 8,
               label_metacluster = 3)
}





