####################################################
#### INSTALLATION OF ALL PACKAGES Below - do not need to repeat if already in library (Try running lines 16-19 first, if some packages missing revert to line 4-11):
####################################################
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("FlowSOM")
install.packages("biocLite")
biocLite(suppressUpdates = TRUE)
biocLite("flowCore", suppressUpdates = TRUE)
install.packages('devtools')
install.packages('Rcpp')

install.packages('biclust')
install.packages('data.table')
install.packages('diptest')
install.packages('evtree')
install.packages('ggdendro')
install.packages("ggfortify")
install.packages('ggplot2')
install.packages('gplots')
install.packages('gdata')
install.packages('ggrepel')
install.packages('ggRandomForests')
install.packages('gridExtra')
install.packages('gtable')
install.packages('gtools')
install.packages('igraph')
install.packages('MASS')
install.packages('packcircles')
install.packages('plyr')
install.packages("randomForestSRC")
install.packages('reshape2')
install.packages('pheatmap')
install.packages('readxl')
install.packages("raster")
install.packages('openxlsx')

install.packages('devtools')
install.packages("ps")
install.packages("backports")
library("backports")

library("ps")
library("devtools")

install_github('tchitchek-lab/SPADEVizR')

source("http://bioconductor.org/biocLite.R")
biocLite(suppressUpdates = TRUE)
biocLite("flowCore", suppressUpdates = TRUE)

install.packages('edgeR')
biocLite("edgeR")
install.packages("bindrcpp")
install.packages("stringi")
install.packages("statmod")
install.packages("tidyselect")
install.packages("tibble")
###################################################
# Library the packages
###################################################
library("devtools")
library("FlowSOM")
library("tidyselect")
library('Rcpp')
library("SPADEVizR")
library(statmod)
library("edgeR")
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(readxl)
library(openxlsx)
library(data.table)
library(ggplot2)
library("tibble")
library(raster)
####################################################
####################################################
source("utils.R") #Sources utils function for phenoviewer_modified
##################################################
# Parallel coordinate plots generated using SPADEvizR - FOR GROUP 1 DATA:
##################################################
### Imports Sheet 4 and Sheet 2 and renames "Abundance" and "Phenotype" respectively, from desired excel file - must change path and excel file name for particular function
PrimaryDirectory <- getwd()
Abundance <- read_excel("./DR3 D52 Grp 1 B Cells 20190103 K=49.xlsx", sheet = "Sheet4")
View(Abundance)
Phenotype <- read_excel("./DR3 D52 Grp 1 B Cells 20190103 K=49.xlsx", sheet = "Sheet2")
View(Phenotype)

### Reformats data for R to run SpadeVizR Script - must change lines 43 and 47 to match size of Abundance and Phenotype Sheets (rows, columns)
cluster.abundances <- as.data.frame(Abundance[1:123,1:17])
rownames(cluster.abundances) <- cluster.abundances[,1]
cluster.abundances <- cluster.abundances[,-1]

cluster.phenotypes <- as.data.frame(Phenotype[1:1968,1:37])
cluster.phenotypes <- cluster.phenotypes[,-3]
results <- importResultsFromTables(cluster.abundances = cluster.abundances, cluster.phenotypes = cluster.phenotypes)


 ### MODIFIED PHENOVIEWER SCRIPT FOR MORE ACCURATE PARALLEL PLOTS ### 
 phenoViewer_modified <- function(Results,
                                  samples        = NULL,
                                  clusters       = NULL,
                                  markers        = NULL,
                                  show.mean      = "both",
                                  show.on_device = TRUE,
                                  sort.markers   = TRUE) {
   
   ### when testing the function, use the parameters inside the function and test line by line of code. Use statement below to test the function above
   # Results=results
   # samples        = NULL
   # clusters       = "Cluster 10"
   # markers        = NULL
   # show.mean      = "only"
   # show.on_device = TRUE
   # sort.markers   = TRUE
   
   if (is.null(Results)) {
     stop("Error in phenoViewer: 'Results' parameter can not be NULL")
   } else if (class(Results)[1] != "Results") {
     stop("Error in phenoViewer: 'Results' parameter must be a 'Results' object")
   }
   
   if(length(Results@marker.names) == 0){
     stop("Error in phenoViewer: 'Results' object must contain phenotypes")
   }
   
   if (is.null(samples)) {
     samples     <- Results@sample.names
     data        <- Results@cluster.phenotypes
     cluster.abundances <- Results@cluster.abundances
   } else if (!all(samples %in% Results@sample.names)) {
     stop("Error in phenoViewer: 'samples' parameter must contains only samples names\n Unknown sample names: ",
          paste(setdiff(unique(samples), Results@sample.names), collapse = " "))
   } else {
     data               <- subset(Results@cluster.phenotypes, sample %in% samples, drop = FALSE)
     cluster.abundances <- Results@cluster.abundances[, samples, drop = FALSE]
   }
   
   data <- stats::na.omit(data)
   
   if (is.null(clusters)) {
     stop("Error in phenoViewer: 'clusters' parameter is required")
   } else if (all(clusters %in% Results@cluster.names)) {
     if (typeof(clusters) != "character") {
       stop("Error in phenoViewer: 'clusters' parameter must be a character vector")
     }
     clusters        <- unique(clusters)
     clusters.select <- data[, "cluster"] %in% clusters
     data            <- data[clusters.select,]
     cluster.abundances     <- cluster.abundances[clusters,]
   } else {
     stop("Error in phenoViewer:\nUnknown clusters : ", paste(setdiff(unique(clusters), Results@cluster.names), collapse = " "))
   }
   
   data <- plyr::ddply(data, c("sample"), function(df) {
     apply(df[, 3:ncol(df)], 2, mean, na.rm = TRUE)
   }) 
   
   if (is.null(markers)) {
     markers <- Results@marker.names
   } else if (all(markers %in% Results@marker.names)) {
     markers <- unique(markers)
     data <- data[, c("sample", markers)]
   } else {
     stop("Error in phenoViewer: Unknown markers :", paste(setdiff(unique(markers), Results@marker.names), collapse = " "))
   }
   
   if (show.mean != "none" && show.mean != "both" && show.mean != "only") {
     stop("Error in phenoViewer: 'show.mean' parameter must contain only one of these : 'none', 'both' or 'only'")
   }
   
   if (!is.logical(show.on_device)) { stop("Error in phenoViewer: 'show.on_device' parameter must be a logical") }
   
   data           <- reshape2::melt(data, id = c("sample"), stringsAsFactors = FALSE)
   colnames(data) <- c("samples", "marker", "value")
   
   names.palette  <- unique(Results@cluster.phenotypes$sample)
   palette        <- ggcolors(length(names.palette))
   names(palette) <- names.palette
   
   assignments <- Results@assignments
   
   if (!is.null(assignments)) {
     
     order       <- unique(assignments$bc)
     assignments <- assignments[samples, , drop = FALSE]
     data$bc <- assignments[data$samples, "bc"]
     order       <- intersect(order, unique(assignments$bc))
     data$bc <- factor(data$bc, levels = order)
     
     names.palette  <- unique(assignments$bc)
     palette        <- ggcolors(length(names.palette))
     names(palette) <- names.palette
     
   } else if (is.element("bc", colnames(assignments))) {
     warning("Warning in phenoViewer: 'assignments' slot do not contain the column 'bc' in the provided 'Results' object. Consequently, the samples names will be used in remplacement")
   } else {
     warning("Warning in phenoViewer: 'assignments' slot in the provided 'Results' object is absent. Consequently, the samples names will be used in remplacement")
   }
   
   if(sort.markers==TRUE){
     clustering.markers  <- Results@clustering.markers
     ordered.markers     <- c(gtools::mixedsort(clustering.markers),gtools::mixedsort(setdiff(Results@marker.names, clustering.markers)))
     bold.markers        <- ifelse(is.element(ordered.markers, clustering.markers), "bold", "plain")
     colored.markers     <- ifelse(is.element(ordered.markers, clustering.markers), "blue", "black")
     data$marker         <- factor(data$marker, levels = ordered.markers, ordered = TRUE)
   }else{
     clustering.markers  <- Results@clustering.markers
     ordered.markers     <- markers
     bold.markers        <- ifelse(is.element(ordered.markers, clustering.markers), "bold", "plain")
     colored.markers     <- ifelse(is.element(ordered.markers, clustering.markers), "blue", "black")
     data$marker         <- factor(data$marker, levels = ordered.markers, ordered = TRUE)
   }
   
   for (i in seq_len(nrow(data))) {
     data[i, "lower.bound"] <- Results@bounds[1, as.character(data[i, "marker"])]
     data[i, "upper.bound"] <- Results@bounds[2, as.character(data[i, "marker"])]
   }
   
   cells.number <- sum(colSums(cluster.abundances))
   
   title    <- paste("Pheno Viewer - cluster: ", paste0(clusters, collapse = ", "), " (", format(cells.number, big.mark = " "), " cells)", sep = "")
   bounds   <- as.numeric(row.names(Results@bounds))
   subtitle <- paste0("Grey ribbon displays from ", (bounds[1] * 100), "% to ", (bounds[2] * 100), "% percentiles of the range expression")
   
   max.value <- -1
   min.value <- -1
   
   max.value <- max(c(data$value, data$upper.bound), na.rm = TRUE)
   min.value <- min(c(data$value, data$lower.bound), na.rm = TRUE)
   
   max.value <-  max.value * (1 + sign(max.value) * 0.1)
   min.value <-  min.value * (1 - sign(min.value) * 0.1)
   
   
   means <- plyr::ddply(data,
                        c("marker"),
                        function(df){mean(df$value, na.rm = TRUE)})
   colnames(means) <- c("marker", "means")
   
   data_means <- data.frame(marker = 0, means= 0, clusters = 0)
   tmp_clusters<- unique(cluster.phenotypes$Cluster) ###### make sure the cluster.phenotypes file column name is "Cluster" and not "cluster"
   for(i in tmp_clusters){
     tmp_data<- Results@cluster.phenotypes
     tmp_clusters.select <- tmp_data[, "cluster"] %in% i
     tmp_data <- tmp_data[tmp_clusters.select,]
     
     tmp_data <- plyr::ddply(tmp_data, c("sample"), function(df) {
       apply(df[, 3:ncol(df)], 2, mean, na.rm = TRUE)
     }) 
     
     tmp_data           <- reshape2::melt(tmp_data, id = c("sample"), stringsAsFactors = FALSE)
     colnames(tmp_data) <- c("samples", "marker", "value")
     
     tmp_means <- plyr::ddply(tmp_data,
                              c("marker"),
                              function(df){mean(df$value, na.rm = TRUE)})
     colnames(tmp_means) <- c("marker", "means")
     tmp_means$clusters = i
     
     data_means = rbind(data_means, tmp_means)
   }
   data_means = data_means[-1, ]
   # data_means$marker = substr(data_means$marker, 2, 100000)
   #data_means = data_means[order(data_means$marker, decreasing = TRUE), ]
   plot <- ggplot2::ggplot(data = data_means) +
     ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), ""))))
   
   plot  <- plot + ggplot2::geom_line(ggplot2::aes_string(x = "marker", y = "means", group = "clusters"),
                                      size = 0.5, #changes size of background lines
                                      alpha = 1,
                                      color = "#CCCCCC")+ 
     ggplot2::scale_y_continuous(limits = c(min.value, max.value), breaks = round(seq(0, max.value, by = 1), 0)) +
     ggplot2::theme_bw()
   
   plot <- plot + ggplot2::geom_line(data  = means,
                                     ggplot2::aes_string(x = "marker", y = "means", group = 1),
                                     #group = 1,
                                     linetype = "solid",
                                     size  = 1,
                                     color = "#FF6666") 
   
   plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, face = bold.markers, color = colored.markers)) +
     ggplot2::theme(legend.text = ggplot2::element_text(size = 6),
                    legend.key  = ggplot2::element_blank(),
                    plot.title  = ggplot2::element_text(hjust=0.5)) +
     ggplot2::xlab("markers") +
     ggplot2::ylab("marker expressions") +
     ggplot2::guides(col = ggplot2::guide_legend(ncol = 1))
   
   
   grid::grid.draw(plot)
   invisible(plot)
 } 

dir.create("Group1_ClusterImages", showWarnings = FALSE)
setwd("Group1_ClusterImages")
for(i in 1:nrow(cluster.abundances)){
  #i=1
  jpeg(paste(rownames(cluster.abundances)[i], ".jpeg", sep = ""),
       width=2000,
       height=1500, 
       res = 300)
  phenoViewer_modified(results, clusters = rownames(cluster.abundances)[i])
  dev.off()
}
setwd(PrimaryDirectory)

####################################
#SCATTER PLOT GENERATOR
###################################
GroupOne_SheetFour <- read_excel("./DR3 D52 Grp 1 B Cells 20190103 K=49.xlsx", sheet = "Sheet4")
write.table(GroupOne_SheetFour, file = "Data for Scatter Plot Group 1.txt", sep = "\t",row.names = FALSE, col.names = TRUE)

#load data
data <- read.table("Data for Scatter Plot Group 1.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data <- as.data.frame(data[1:123,1:17])
data$Cluster <- gsub(" ", "_", data$Cluster, fixed = TRUE)
rownames(data) <- data$Cluster
data <- data[,-1]

sum_counts_sample <- colSums(data)

for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
    data[i,j] = data[i,j]/sum_counts_sample[j]*100
  }
}

#transpose the data for ploting
data <- t(data)
data <- as.data.frame(data)

#group assignment
group_data <- read.table("group assignment for group 1.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
group_data$sample <- trim(group_data$sample)
group_data$sample = gsub(" ", ".", group_data$sample, fixed = TRUE)
data$group <- group_data$group#[match(rownames(data), group_data$sample)]
data <- data[, c(ncol(data), 1:(ncol(data)-1))]

dir.create("Group1_Scatterplots", showWarnings = FALSE)
setwd("Group1_Scatterplots")
x_order = factor(data$group, levels=c("Naive Unchallenged", "Naive Challenged", "Vax Unchallenged","Vax Challenged"), ordered=TRUE)
for(i in 2:ncol(data)){
  scatter_plot <- 
    ggplot(data, aes_string(x = x_order, fill = "group", y = colnames(data)[i]))+ 
    geom_dotplot(binaxis = "y", stackdir = "centerwhole") +
    stat_summary(fun.y = "median", size=0.5, geom = 'line', aes(group=1))+
    stat_summary(
      fun.ymin = function(z) { quantile(z,0.25) },
      fun.ymax = function(z) { quantile(z,0.75) },
      fun.y = median,
      width = 0.2,
      geom = "errorbar") + 
    theme(axis.text.x = element_text(size = 25, face = "bold", vjust = 1.0, hjust = 1.0, angle = 45)) +
    theme(axis.text.y = element_text(size = 25, face = "bold", vjust = 0.5, hjust = 0.5, angle = 0)) +
    theme(legend.position = "none")
  ggsave(scatter_plot,
         width = 20,
         height = 15,
         dpi = 300,
         filename = paste(colnames(data)[i], ".jpeg", sep = ""))
}
setwd(PrimaryDirectory)

##################################################
# Parallel co-ordinate plots generated using SPADEvizR  - FOR GROUP 2 DATA:
##################################################
### Imports Sheet 4 and Sheet 2 and renames "Abundance" and "Phenotype" respectively, from desired excel file - must change path and excel file name for particular function
Abundance <- read_excel("./DR3 D52 Grp 2 B Cells 20190103 K=41.xlsx", sheet = "Sheet4")
View(Abundance)
Phenotype <- read_excel("./DR3 D52 Grp 2 B Cells 20190103 K=41.xlsx", sheet = "Sheet2")
View(Phenotype)

### Reformats data for R to run SpadeVizR Script - must change lines 85 and 89 to match size of Abundance and Phenotype Sheets (rows, columns)
cluster.abundances <- as.data.frame(Abundance[1:117,1:16])
rownames(cluster.abundances) <- cluster.abundances[,1]
cluster.abundances <- cluster.abundances[,-1]

cluster.phenotypes <- as.data.frame(Phenotype[1:1755,1:37])
cluster.phenotypes <- cluster.phenotypes[,-3]
results <- importResultsFromTables(cluster.abundances = cluster.abundances, cluster.phenotypes = cluster.phenotypes)


dir.create("Group2_ClusterImages", showWarnings = FALSE)
setwd("Group2_ClusterImages")
for(i in 1:nrow(cluster.abundances)){
  jpeg(paste(rownames(cluster.abundances)[i], ".jpeg", sep = ""),
       width=2000,
       height=1500, 
       res = 300)
  phenoViewer_modified(results, clusters = rownames(cluster.abundances)[i])
  dev.off()
}
setwd(PrimaryDirectory)

####################################
#SCATTER PLOT GENERATOR
###################################
GroupTwo_SheetFour <- read_excel("./DR3 D52 Grp 2 B Cells 20190103 K=41.xlsx", sheet = "Sheet4")
write.table(GroupTwo_SheetFour, file = "Data for Scatter Plot Group 2.txt", sep = "\t",row.names = FALSE, col.names = TRUE)

#load data
data <- read.table("Data for Scatter Plot Group 2.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data <- as.data.frame(data[1:117,1:16])
data$Cluster <- gsub(" ", "_", data$Cluster, fixed = TRUE)
rownames(data) <- data$Cluster
data <- data[,-1]

sum_counts_sample <- colSums(data)

for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
    data[i,j] = data[i,j]/sum_counts_sample[j]*100
  }
}

#transpose the data for ploting
data <- t(data)
data <- as.data.frame(data)

#group assignment
group_data <- read.table("group assignment for group 2.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
group_data$sample <- trim(group_data$sample)
group_data$sample = gsub(" ", ".", group_data$sample, fixed = TRUE)
data$group <- group_data$group#[match(rownames(data), group_data$sample)]
data <- data[, c(ncol(data), 1:(ncol(data)-1))]

dir.create("Group2_Scatterplots", showWarnings = FALSE)
setwd("Group2_Scatterplots")
x_order = factor(data$group, levels=c("Naive Unchallenged", "Naive Challenged", "Vax Unchallenged","Vax Challenged"), ordered=TRUE)
for(i in 2:ncol(data)){
  scatter_plot <- 
    ggplot(data, aes_string(x = x_order, fill = "group", y = colnames(data)[i]))+ 
    geom_dotplot(binaxis = "y", stackdir = "centerwhole") +
    stat_summary(fun.y = "median", size=0.5, geom = 'line', aes(group=1))+
    stat_summary(
      fun.ymin = function(z) { quantile(z,0.25) },
      fun.ymax = function(z) { quantile(z,0.75) },
      fun.y = median,
      width = 0.2,
      geom = "errorbar") + 
    theme(axis.text.x = element_text(size = 25, face = "bold", vjust = 1.0, hjust = 1.0, angle = 45)) +
    theme(axis.text.y = element_text(size = 20, face = "bold", vjust = 0.5, hjust = 0.5, angle = 0)) +
    theme(legend.position = "none")
  ggsave(scatter_plot,
         width = 20,
         height = 15,
         dpi = 300,
         filename = paste(colnames(data)[i], ".jpeg", sep = ""))
}
setwd(PrimaryDirectory)

#########################################################################################################################################################
#########################################################################################################################################################
#### Next four lines of code generate .txt files from sheet one of group 1 and 2 excel sheets to be used for pearson's correlation
### Creates a .txt file of sheet one from Group 1 excel file containing all vortex data
GroupOne_SheetOne <- read_excel("./DR3 D52 Grp 1 B Cells 20190103 K=49.xlsx", sheet = "Sheet1")
write.table(GroupOne_SheetOne, file = "DR3 D52 Grp 1 B Cells 20190103 K=49.txt", sep = "\t",row.names = FALSE, col.names = TRUE)

### Creates a .txt file of sheet one from Group 2 excel file containing all vortex data
GroupTwo_SheetOne <- read_excel("./DR3 D52 Grp 2 B Cells 20190103 K=41.xlsx", sheet = "Sheet1")
write.table(GroupTwo_SheetOne, file = "DR3 D52 Grp 2 B Cells 20190103 K=41.txt", sep = "\t",row.names = FALSE, col.names = TRUE)


###################################################
# Generates a list of matching clusters from group 1 and 2 based on pearson's correlation and count
###################################################
rescale_to_0_1 <- function(experiment_name, experiment_file, rescale = TRUE){

  #read the file
  raw_table = read.delim(experiment_file, sep = "\t", stringsAsFactors = FALSE)
  
  #modify the column name
  colnames(raw_table) = gsub("X.", "", colnames(raw_table), fixed = TRUE)
  colnames(raw_table) = gsub("X", "", colnames(raw_table), fixed = TRUE)
  
  #modify the cluster name
  raw_table$Cluster = gsub("Cluster", "", raw_table$Cluster, fixed = TRUE)
  raw_table$Cluster = gsub(" ", "", raw_table$Cluster, fixed = TRUE)
  raw_table$Cluster = paste("Cluster_", raw_table$Cluster, sep = "")
  
  #sorting the dataset for better view
  raw_table = raw_table[order(raw_table$Cluster, raw_table$Term, decreasing = FALSE),]
  
  #obtain the samples name
  samples = unique(raw_table$Term)
  
  #obtain the amount of samples
  nSample= length(samples)
  
  #obtain the cluster name
  clusters = unique(raw_table$Cluster)
  
  #obtain the amount of clusters
  nCluster = length(clusters)
  
  #create a blank table with labels
  mean_marker_total_cells = data.frame(tmp_name = 0)
  
  for(i in 1:(ncol(raw_table)-1)){
    mean_marker_total_cells = cbind(mean_marker_total_cells, 0)
  }
  
  colnames(mean_marker_total_cells) = colnames(raw_table)
  mean_marker_total_cells = mean_marker_total_cells[, colnames(mean_marker_total_cells)!= "Term"]
  
  # creat a array for storing the number of total numbers of cluster for samples
  cluster_count = rep(0,nCluster)
  
  #calculate and store the total numbers of cluster for samples
  j = 1
  k = 1
  for(i in 1:nrow(raw_table)){
    if(i == nrow(raw_table)){
      cluster_count[j] = cluster_count[j] + 1
      for(n in 1:ncol(mean_marker_total_cells)){
        if(n == 1){
          mean_marker_total_cells[k,n] = raw_table$Cluster[i]
        }
        if(n == 2){
          mean_marker_total_cells[k,n] = sum(raw_table[(i-cluster_count[j]+1):i,n+1])
        }
        if(n > 2){
          mean_marker_total_cells[k,n] = mean(raw_table[(i-cluster_count[j]+1):i,n+1])
        }
      }
      break()
    }
    
    if(raw_table$Cluster[i] == raw_table$Cluster[i+1]){
      cluster_count[j] = cluster_count[j] + 1
    }else{
      cluster_count[j] = cluster_count[j] + 1
      for(n in 1:ncol(mean_marker_total_cells)){
        if(n == 1){
          mean_marker_total_cells[k,n] = raw_table$Cluster[i]
        }
        if(n == 2){
          mean_marker_total_cells[k,n] = sum(raw_table[(i-cluster_count[j]+1):i,n+1])
        }
        if(n > 2){
          mean_marker_total_cells[k,n] = mean(raw_table[(i-cluster_count[j]+1):i,n+1])
        }
      }
      mean_marker_total_cells = rbind(mean_marker_total_cells, 0)
      j = j + 1
      k = k + 1
    }
  }
  
  tmp_rescale <- function(x) (x-min(x))/(max(x) - min(x))
  
  tmp_mean_marker_total_cells = mean_marker_total_cells
  
  tmp_mean_marker_total_cells$Cluster = paste(experiment_name, 
                                              "_",
                                              tmp_mean_marker_total_cells$Cluster, 
                                              sep = "")
  
  rownames(tmp_mean_marker_total_cells) = tmp_mean_marker_total_cells[,1]
  tmp_mean_marker_total_cells = tmp_mean_marker_total_cells[,-1]
  
  tmp_mean_marker_total_cells$Count1 = tmp_mean_marker_total_cells$Count
  
  if(rescale==TRUE){
    for(i in 1:(ncol(tmp_mean_marker_total_cells)-1)){
      tmp_mean_marker_total_cells[,i] = tmp_rescale(tmp_mean_marker_total_cells[,i])
    }
  }
  
  return(tmp_mean_marker_total_cells)
} 

### Change experiment_file names to match reformatted excel sheets used for SpadevizR
## All files must be in correct folder in the working path in order to run code!
experiment1 = rescale_to_0_1(experiment_name = "Grp1",
                             experiment_file = "DR3 D52 Grp 1 B Cells 20190103 K=49.txt",
                             rescale = TRUE)
#includes count - can # if dont want to rank based on count
#experiment1 = experiment1[, colnames(experiment1) != "Count"] 

experiment2 = rescale_to_0_1(experiment_name = "Grp2",
                             experiment_file = "DR3 D52 Grp 2 B Cells 20190103 K=41.txt",
                             rescale = TRUE)
#includes count - can # if dont want to rank based on count
#experiment2 = experiment2[, colnames(experiment2) != "Count"] 

experiment1_1 = experiment1
experiment1 = experiment1[, colnames(experiment1) != "Count1"]

experiment2_1 = experiment2
experiment2 = experiment2[, colnames(experiment2) != "Count1"]

#create a blank table to store the pearson correlation results
experiment1_experiment2_Pearson_correlation<-data.frame(experiment1_cluster = 0, experiment2_cluster = 0, experiment1_count = 0, experiment2_count = 0)

#perform pairwise pearson correlation between experiment1 and experiment2
t=1
for(i in 1:nrow(experiment1)){
  for(j in 1:nrow(experiment2)){
    experiment1_experiment2_Pearson_correlation$experiment1_cluster[t]<-rownames(experiment1)[i]
    experiment1_experiment2_Pearson_correlation$experiment2_cluster[t]<-rownames(experiment2)[j]
    
    pearson_statictis<-cor.test(as.numeric(experiment1[i,]),as.numeric(experiment2[j,]),method = "pearson")
    
    experiment1_experiment2_Pearson_correlation$cor[t]<-pearson_statictis$estimate
    experiment1_experiment2_Pearson_correlation$p.value[t]<-pearson_statictis$p.value
    
    experiment1_experiment2_Pearson_correlation$experiment1_count[t]<-experiment1_1$Count1[i]
    experiment1_experiment2_Pearson_correlation$experiment2_count[t]<-experiment2_1$Count1[j]
    
    t<-t+1
    experiment1_experiment2_Pearson_correlation<-rbind(experiment1_experiment2_Pearson_correlation, 0)
  }
}
experiment1_experiment2_Pearson_correlation = experiment1_experiment2_Pearson_correlation[, c(1,2,5,6,3,4)]

#Sorting the data for better view
experiment1_experiment2_Pearson_correlation = experiment1_experiment2_Pearson_correlation[
  order(experiment1_experiment2_Pearson_correlation$cor, decreasing = TRUE),]

#Take a look at the results
View(experiment1_experiment2_Pearson_correlation)

# create a CSV file storing the pearson correlation data
write.csv(experiment1_experiment2_Pearson_correlation, "DR3 D52 B Cells Pearsons Coefficient.csv", row.names = FALSE)

####################################################################################################################################################################################
####################################################################################################################################################################################

setwd(PrimaryDirectory)

##################################################
## Matching Script - combines cluster data from group 1 and 2 and reformats to create grouped_file containing all cluster information for newly named matched clusters (n=4 -> n=8)
##################################################

# read data
Grp1_file <- "DR3 D52 Grp 1 B Cells 20190103 K=49.txt"
Grp1_data <- read.table(Grp1_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Grp1_data$Cluster = paste("Grp1_Cluster_", Grp1_data$Cluster, sep ="")

Grp2_file <- "DR3 D52 Grp 2 B Cells 20190103 K=41.txt"
Grp2_data <- read.table(Grp2_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Grp2_data$Cluster = paste("Grp2_Cluster_", Grp2_data$Cluster, sep ="")

for(i in 4:37){
  #i = 4
  Grp1_min = min(Grp1_data[,i])
  Grp2_min = min(Grp2_data[,i])
  if(Grp1_min > Grp2_min){
    correlation = Grp1_min - Grp2_min
    Grp1_data[,i] = Grp1_data[,i] - correlation
  }
  if(Grp2_min > Grp1_min){
    correlation = Grp2_min - Grp1_min
    Grp2_data[,i] = Grp2_data[,i] - correlation
  }
}

grouped_file = rbind(Grp1_data, Grp2_data)
### MAKE SURE THE COLUMN V ("CXCR4" has the X in it. Often group 2 sheet will read "CCR4")

# Change matching_file_name to name of .txt file that contains Group 1 Clusters and their new_name in addition to matching Group 2 Clusters and their new_name  
matching_file_name = "Matched Clusters DR3 Day 52 B Cell.txt"

#grouped_file = read.table(grouped_file_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
matching_file = read.table(matching_file_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

matching_file = cbind(matching_file, 0)
colnames(matching_file)[3] = "Cluster_new"


for(i in 1:(nrow(matching_file)/2)){
  matching_file$Cluster_new[i] = paste("Grp1_Cluster_", matching_file$Cluster[i], sep = "")
}

for(i in ((nrow(matching_file)/2)+1):nrow(matching_file)){
  matching_file$Cluster_new[i] = paste("Grp2_Cluster_", matching_file$Cluster[i], sep = "")
}

matching_file = matching_file[,c(3,1,2)]
matching_file$Cluster = NULL
colnames(matching_file)[1] = "Cluster"

#match the new name
grouped_file$new_cluster_name = matching_file$new_name[match(grouped_file$Cluster, matching_file$Cluster)]

#resort the data
grouped_file = grouped_file[, c(1, 2, ncol(grouped_file), 3:(ncol(grouped_file)-1))]

#Replace NA in new_cluster_name with original cluster number
for (i in 1: nrow(grouped_file)){
  if(is.na(grouped_file$new_cluster_name[i]) == TRUE){
    grouped_file$new_cluster_name[i] = substr(grouped_file$Cluster[i], 1, 20)
  }
}

# #delete the unmatched clusters
# grouped_file = grouped_file[!is.na(grouped_file$new_cluster_name), ]

#sorting the data for better view
grouped_file = grouped_file[order(grouped_file$new_cluster_name),]


grouped_file = cbind(grouped_file, 0)
colnames(grouped_file)[39] = "Cluster_new"
grouped_file = grouped_file[,c(39,1:ncol(grouped_file))]
grouped_file$Cluster = NULL
colnames(grouped_file)[1] = "Cluster"
grouped_file$Cluster = grouped_file$new_cluster_name
grouped_file$new_cluster_name = NULL
grouped_file$Cluster_new.1 = NULL



write.xlsx(grouped_file, "Phenotype DR3 Day52 B Cells ALL BACKGROUND.xlsx", row.names=FALSE)

##############################################
# Generates an "Abundance" sheet to use for analysis using the frequency of the cluster in the mouse 
##############################################
#-------recale function-------#
tmp_percent <- function(x) (x/sum(x))*100
#-------recale function-------#


# Change grouped_file_name to name of .txt file that contains ALL DATA from Group 1 and 2

Grp1_file = "Data for Scatter Plot Group 1.txt"
Grp2_file = "Data for Scatter Plot Group 2.txt"

Grp1_data = read.table(Grp1_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Grp1_data <- as.data.frame(Grp1_data[1:123,1:17])
Grp1_data$Cluster = paste("Grp1_", Grp1_data$Cluster, sep ="")
Grp1_data$Cluster = gsub(" ", "_", Grp1_data$Cluster)
#rescale the markers expression
for(i in 2:ncol(Grp1_data)){
  Grp1_data[,i] = tmp_percent(Grp1_data[,i])
}


Grp2_data = read.table(Grp2_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Grp2_data <- as.data.frame(Grp2_data[1:117,1:16])
Grp2_data$Cluster = paste("Grp2_", Grp2_data$Cluster, sep ="")
Grp2_data$Cluster = gsub(" ", "_", Grp2_data$Cluster)
#rescale the markers expression
for(i in 2:ncol(Grp2_data)){
  Grp2_data[,i] = tmp_percent(Grp2_data[,i])
}


matching_file_name = "Matched Clusters DR3 Day 52 B Cell.txt"

#grouped_file = read.table(grouped_file_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
matching_file = read.table(matching_file_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

matching_file = cbind(matching_file, 0)
colnames(matching_file)[3] = "Cluster_new"


for(i in 1:(nrow(matching_file)/2)){
  matching_file$Cluster_new[i] = paste("Grp1_Cluster_", matching_file$Cluster[i], sep = "")
}

for(i in ((nrow(matching_file)/2)+1):nrow(matching_file)){
  matching_file$Cluster_new[i] = paste("Grp2_Cluster_", matching_file$Cluster[i], sep = "")
}

matching_file = matching_file[,c(3,1,2)]
matching_file$Cluster = NULL
colnames(matching_file)[1] = "Cluster"


#create a blank table to store data
grouped_data = data.frame(tmp_name = 0)
tmp_data = matching_file[matching_file$new_name == 1,]

Grp1_tmp_data = Grp1_data[Grp1_data$Cluster %in% tmp_data$Cluster, ]
rownames(Grp1_tmp_data) = Grp1_tmp_data$Cluster
Grp1_tmp_data$Cluster = NULL

Grp2_tmp_data = Grp2_data[Grp2_data$Cluster %in% tmp_data$Cluster, ]
rownames(Grp2_tmp_data) = Grp2_tmp_data$Cluster
Grp2_tmp_data$Cluster = NULL

tmp_combined_data_1 = cbind(Grp1_tmp_data, Grp2_tmp_data)
#tmp_combined_data_1 = rbind(Grp1_tmp_data, Grp2_tmp_data)
#rownames(tmp_combined_data_1) = paste("Cluster_", 1, sep = "")

grouped_data = tmp_combined_data_1

for (i in unique(matching_file$new_name)){
  #i=1
  
  tmp_data = matching_file[matching_file$new_name == i,]
  
  Grp1_tmp_data = Grp1_data[Grp1_data$Cluster %in% tmp_data$Cluster, ]
  rownames(Grp1_tmp_data) = Grp1_tmp_data$Cluster
  Grp1_tmp_data$Cluster = NULL
  
  Grp2_tmp_data = Grp2_data[Grp2_data$Cluster %in% tmp_data$Cluster, ]
  rownames(Grp2_tmp_data) = Grp2_tmp_data$Cluster
  Grp2_tmp_data$Cluster = NULL  
  
  tmp_combined_data = cbind(Grp1_tmp_data, Grp2_tmp_data)
  rownames(tmp_combined_data) = paste("Cluster_", i, sep = "")
  
  grouped_data = rbind(grouped_data, tmp_combined_data)
  
}

grouped_data = grouped_data[-1, ]

#grouped_file = grouped_file[order(grouped_file$new_cluster_name),]

write.xlsx(grouped_data, "Abundance DR3 Day 52 B Cell Data.xlsx", row.names=TRUE)

#####################################################################################################
#MODIFIED PCP GENERATOR FOR PHENOTYPE
#####################################################################################################
setwd(PrimaryDirectory)
### Import "Phenotype", from desired excel file - must change path and excel file name for particular function
##### MAKE SURE YOU HAVE SAVED YOUR GENERATED PHENOTYPE SHEET AS AN .XLSX File
Phenotype <- read_excel("./Phenotype DR3 Day52 B Cells ALL BACKGROUND.xlsx", sheet = "Sheet 1")
View(Phenotype)

#Must change parameters of phenotype sheet according to file size
cluster.phenotypes <- as.data.frame(Phenotype[1:3723,1:37])
cluster.phenotypes <- cluster.phenotypes[,-3]

phenoViewer_modified_v2 <-function(  cluster.phenotypes,
                                     samples        = NULL,
                                     clusters,
                                     markers        = NULL,
                                     show.mean      = "only",
                                     show.on_device = TRUE,
                                     sort.markers   = TRUE){
  
  if (is.null(samples)) {
    samples <- unique(cluster.phenotypes$Term)
    data        <- cluster.phenotypes
  } else if (!all(samples %in% Results@sample.names)) {
    stop("Error in phenoViewer: 'samples' parameter must contains only samples names\n Unknown sample names: ",
         paste(setdiff(unique(samples), Results@sample.names), collapse = " "))
  } else {
    data               <- subset(Results@cluster.phenotypes, sample %in% samples, drop = FALSE)
    cluster.abundances <- Results@cluster.abundances[, samples, drop = FALSE]
  }
  
  data <- stats::na.omit(data)
  
  clusters        <- unique(clusters)
  clusters.select <- data[, "Cluster"] %in% clusters
  data            <- data[clusters.select,]
  
  data <- plyr::ddply(data, c("Term"), function(df) {
    apply(df[, 3:ncol(df)], 2, mean, na.rm = TRUE)
  }) 
  
  data           <- reshape2::melt(data, id = c("Term"), stringsAsFactors = FALSE)
  colnames(data) <- c("samples", "marker", "value")
  
  title    <- paste("Cluster_", clusters, sep = "")
  
  max.value <- -1
  min.value <- -1
  
  max.value <- max(c(data$value, data$upper.bound), na.rm = TRUE)
  min.value <- min(c(data$value, data$lower.bound), na.rm = TRUE)
  
  max.value <-  max.value * (1 + sign(max.value) * 0.1)
  min.value <-  min.value * (1 - sign(min.value) * 0.1)
  
  means <- plyr::ddply(data,
                       c("marker"),
                       function(df){mean(df$value, na.rm = TRUE)})
  colnames(means) <- c("marker", "means")
  
  data_means <- data.frame(marker = 0, means= 0, clusters = 0)
  tmp_clusters<- unique(cluster.phenotypes$Cluster) ###### make sure the cluster.phenotypes file column name is "Cluster" and not "cluster"
  for(i in tmp_clusters){
    tmp_data<- cluster.phenotypes
    tmp_clusters.select <- tmp_data[, "Cluster"] %in% i
    tmp_data <- tmp_data[tmp_clusters.select,]
    
    tmp_data <- plyr::ddply(tmp_data, c("Term"), function(df) {
      apply(df[, 3:ncol(df)], 2, mean, na.rm = TRUE)
    }) 
    
    tmp_data           <- reshape2::melt(tmp_data, id = c("Term"), stringsAsFactors = FALSE)
    colnames(tmp_data) <- c("samples", "marker", "value")
    
    tmp_means <- plyr::ddply(tmp_data,
                             c("marker"),
                             function(df){mean(df$value, na.rm = TRUE)})
    colnames(tmp_means) <- c("marker", "means")
    tmp_means$clusters = i
    
    data_means = rbind(data_means, tmp_means)
  }
  data_means = data_means[-1, ]
  
  rescale_data_means = data_means
  rescale_means = data_means[data_means$clusters == clusters,]
  
  plot <- ggplot2::ggplot(data = rescale_data_means) +
    ggplot2::ggtitle(bquote(atop(.(title))))
  
  plot  <- plot + ggplot2::geom_line(ggplot2::aes_string(x = "marker", y = "means", group = "clusters"),
                                     size = 0.4,
                                     alpha = 1,
                                     color = "#CCCCCC")+ 
    ggplot2::scale_y_continuous(limits = c(min(data_means$means), max(data_means$means)), breaks = round(seq(0, max(data_means$means), by = 1), 0)) +
    ggplot2::theme_bw()
  
  plot <- plot + ggplot2::geom_line(data  = rescale_means,
                                    ggplot2::aes_string(x = "marker", y = "means", group = 1),
                                    linetype = "solid",
                                    size  = 1,
                                    color = "#FF6666") 
  
  plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold")) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 6),
                   legend.key  = ggplot2::element_blank(),
                   plot.title  = ggplot2::element_text(hjust=0.5)) +
    ggplot2::xlab("markers") +
    ggplot2::ylab("marker expressions") +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1))
  
  
  grid::grid.draw(plot)
  invisible(plot)
}

# Should you only want to see one cluster image
# phenoViewer_modified_v2(cluster.phenotypes = cluster.phenotypes,
#                         clusters = "2887")

dir.create("Grouped_ClusterImages", showWarnings = FALSE)
setwd("Grouped_ClusterImages")

a = cluster.phenotypes[which(cluster.phenotypes$Cluster %in% c(1:31)), ]

for (i in unique(a$Cluster)){
  jpeg(paste("Cluster_", i, ".jpeg"),
       width=2000,
       height=1500, 
       res = 300)
  
  phenoViewer_modified_v2(cluster.phenotypes = cluster.phenotypes,
                          clusters = i)
  dev.off()
}

setwd(PrimaryDirectory)
##################################################
# GENERATES PHENOTYPE SHEET FOR GROUPED CLUSTERS TO BE USED FOR SPADEVIZR ANALYSIS
##################################################

# read data
Grp1_file <- "DR3 D52 Grp 1 B Cells 20190103 K=49.txt"
Grp1_data <- read.table(Grp1_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Grp1_data$Cluster = paste("Grp1_Cluster_", Grp1_data$Cluster, sep ="")

Grp2_file <- "DR3 D52 Grp 2 B Cells 20190103 K=41.txt"
Grp2_data <- read.table(Grp2_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Grp2_data$Cluster = paste("Grp2_Cluster_", Grp2_data$Cluster, sep ="")

for(i in 4:37){
  #i = 4
  Grp1_min = min(Grp1_data[,i])
  Grp2_min = min(Grp2_data[,i])
  if(Grp1_min > Grp2_min){
    correlation = Grp1_min - Grp2_min
    Grp1_data[,i] = Grp1_data[,i] - correlation
  }
  if(Grp2_min > Grp1_min){
    correlation = Grp2_min - Grp1_min
    Grp2_data[,i] = Grp2_data[,i] - correlation
  }
}

grouped_file = rbind(Grp1_data, Grp2_data)

# Change matching_file_name to name of .txt file that contains Group 1 Clusters and their new_name in addition to matching Group 2 Clusters and their new_name  
matching_file_name = "Matched Clusters DR3 Day 52 B Cell.txt"

#grouped_file = read.table(grouped_file_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
matching_file = read.table(matching_file_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

matching_file = cbind(matching_file, 0)
colnames(matching_file)[3] = "Cluster_new"


for(i in 1:(nrow(matching_file)/2)){
  matching_file$Cluster_new[i] = paste("Grp1_Cluster_", matching_file$Cluster[i], sep = "")
}

for(i in ((nrow(matching_file)/2)+1):nrow(matching_file)){
  matching_file$Cluster_new[i] = paste("Grp2_Cluster_", matching_file$Cluster[i], sep = "")
}

matching_file = matching_file[,c(3,1,2)]
matching_file$Cluster = NULL
colnames(matching_file)[1] = "Cluster"


#match the new name
grouped_file$new_cluster_name = matching_file$new_name[match(grouped_file$Cluster, matching_file$Cluster)]

#resort the data
grouped_file = grouped_file[, c(1, 2, ncol(grouped_file), 3:(ncol(grouped_file)-1))]

#Replace NA in new_cluster_name with original cluster number
for (i in 1: nrow(grouped_file)){
  if(is.na(grouped_file$new_cluster_name[i]) == TRUE){
    grouped_file$new_cluster_name[i] = substr(grouped_file$Cluster[i], 1, 9)
  }
}

#delete the unmatched clusters
grouped_file = grouped_file[grouped_file$new_cluster_name != "Grp1_Clus", ]
grouped_file = grouped_file[grouped_file$new_cluster_name != "Grp2_Clus", ]

grouped_file = cbind(grouped_file, 0)
colnames(grouped_file)[39] = "Cluster_new"
grouped_file = grouped_file[,c(39,1:ncol(grouped_file))]
grouped_file$Cluster = NULL
colnames(grouped_file)[1] = "Cluster"
grouped_file$Cluster = grouped_file$new_cluster_name
grouped_file$new_cluster_name = NULL
grouped_file$Cluster_new.1 = NULL



#sorting the data for better view
#grouped_file = grouped_file[order(grouped_file$new_cluster_name),]

write.xlsx(grouped_file, "Phenotype DR3 Day52 B Cells SPADEVIZR.xlsx", row.names=FALSE)

# FILE GENERATED ABOVE SERVES AS PHENOTYPE SHEET FOR SPADEVIZR ANALYSIS 


##################################################
# SPADEVIZR ANALYSIS - FOR COMBINED GROUP DATA:
##################################################
### Imports Sheet 4 and Sheet 2 and renames "Abundance" and "Phenotype SpadeVizR" respectively, from desired excel file - must change path and excel file name for particular function
Abundance <- read_excel("./Abundance DR3 Day 52 B Cell Data.xlsx", sheet = "Sheet 1")
View(Abundance)
Phenotype <- read_excel("./Phenotype DR3 Day52 B Cells SPADEVIZR.xlsx", sheet = "Sheet 1")
View(Phenotype)

### Reformats data for R to run SpadeVizR Script - must change lines 334 and 338 to match size of Abundance and Phenotype Sheets (rows, columns)
cluster.abundances <- as.data.frame(Abundance[1:19,1:32])
rownames(cluster.abundances) <- cluster.abundances[,1]
cluster.abundances <- cluster.abundances[,-1]

cluster.phenotypes <- as.data.frame(Phenotype[1:589,1:37])
cluster.phenotypes <- cluster.phenotypes[,-3]
cluster.phenotypes$Cluster = paste("Cluster_", cluster.phenotypes$Cluster, sep ="")
results <- importResultsFromTables(cluster.abundances = cluster.abundances, cluster.phenotypes = cluster.phenotypes)

### Edit file names for each group based on experiment layout (can copy and paste group names from console window below to assure names are correct)
Control <- c("c05_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c06_GRP2_DR3_D52.spleen_B.Cells._37p", "c06_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c07_GRP2_DR3_D52.spleen_B.Cells._37p", "c07_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c08_GRP2_DR3_D52.spleen_B.Cells._37p", "c08_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c09_GRP2_DR3_D52.spleen_B.Cells._37p")
Naive_Chal <- c("c01_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c02_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c03_GRP2_DR3_D52.spleen_B.Cells._37p", "c03_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c04_GRP2_DR3_D52.spleen_B.Cells._37p", "c04_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c05_GRP2_DR3_D52.spleen_B.Cells._37p")
Vax_Chal <- c("c09_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c10_GRP2_DR3_D52.spleen_B.Cells._37p", "c10_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c11_GRP2_DR3_D52.spleen_B.Cells._37p", "c11_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c12_GRP2_DR3_D52.spleen_B.Cells._37p", "c12_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c13_GRP2_DR3_D52.spleen_B.Cells._37p")
Vax_Unchal <- c("c13_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c14_GRP2_DR3_D52.spleen_B.Cells._37p", "c14_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c15_GRP2_DR3_D52.spleen_B.Cells._37p", "c15_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c16_GRP2_DR3_D52.spleen_B.Cells._37p", "c16_Grp1_DR3_Day52_Spleene_B.Cells._37p", "c17_GRP2_DR3_D52.spleen_B.Cells._37p")

### Generates Volcano plots for all conditions selected
## If want to change p-value to 0.01, change "th.pvalue = 0.01"  
# To run an unpaired T-test, method.paired = FALSE. To run a paired T-test use, method.paired = TRUE

### Generates CSV files for all p values for all clusters and saves them in a folder in your working directory
dir.create("SpadevizR Analysis and Volcano Plots", showWarnings = FALSE)
setwd("SpadevizR Analysis and Volcano Plots")

resultsDAC_CvNC <- identifyDAC(results, condition1 = Control, condition2 = Naive_Chal, th.pvalue = 0.05, th.fc = 1, method.paired = FALSE, use.percentages = FALSE)
resultsDAC_CvNC@results
View(resultsDAC_CvNC@results)
write.csv(resultsDAC_CvNC@results, "Control_v_Naive_Chal_DAC_p_values.csv", row.names = FALSE)

tiff("Control vs Naive_chal.tiff",
     width=2000,
     height=1500, 
     res = 300)
SPADEVizR::plot(resultsDAC_CvNC)
dev.off()

resultsDAC_CvVU <- identifyDAC(results, condition1 = Control, condition2 = Vax_Unchal, th.pvalue = 0.05, th.fc = 1, method.paired = FALSE, use.percentages = FALSE)
resultsDAC_CvVU@results
View(resultsDAC_CvVU@results)
write.csv(resultsDAC_CvVU@results, "Control_v_Vax_Unchal_DAC_p_values.csv", row.names = FALSE)

tiff("Control vs Vax_Unchal.tiff",
     width=2000,
     height=1500, 
     res = 300)
SPADEVizR::plot(resultsDAC_CvVU)
dev.off()

resultsDAC_NCvVC <- identifyDAC(results, condition1 = Naive_Chal, condition2 = Vax_Chal, th.pvalue = 0.05, th.fc = 1, method.paired = FALSE, use.percentages = FALSE)
resultsDAC_NCvVC@results
View(resultsDAC_NCvVC@results)
write.csv(resultsDAC_NCvVC@results, "Naive_Chal_v_Vax_Chal_DAC_p_values.csv", row.names = FALSE)

tiff("Naive_Chal vs Vax_Chal.tiff",
     width=2000,
     height=1500, 
     res = 300)
SPADEVizR::plot(resultsDAC_NCvVC)
dev.off()

resultsDAC_VUvVC <- identifyDAC(results, condition1 = Vax_Unchal, condition2 = Vax_Chal, th.pvalue = 0.05, th.fc = 1, method.paired = FALSE, use.percentages = FALSE)
resultsDAC_VUvVC@results
View(resultsDAC_VUvVC@results)
write.csv(resultsDAC_VUvVC@results, "Vax_Unchal_v_Vax_Chal_DAC_p_values.csv", row.names = FALSE)

tiff("Vax_Unchal vs Vax_Chal.tiff",
     width=2000,
     height=1500, 
     res = 300)
SPADEVizR::plot(resultsDAC_VUvVC)
dev.off()


resultsDAC_CvVC <- identifyDAC(results, condition1 = Control, condition2 = Vax_Chal, th.pvalue = 0.05, th.fc = 1, method.paired = FALSE, use.percentages = FALSE)
resultsDAC_CvVC@results
View(resultsDAC_CvVC@results)
write.csv(resultsDAC_CvVC@results, "Control_v_Vax_Chal_DAC_p_values.csv", row.names = FALSE)

tiff("Control vs Vax_Chal.tiff",
     width=2000,
     height=1500, 
     res = 300)
SPADEVizR::plot(resultsDAC_CvVC)
dev.off()

resultsDAC_NCvVU <- identifyDAC(results, condition1 = Naive_Chal, condition2 = Vax_Unchal, th.pvalue = 0.05, th.fc = 1, method.paired = FALSE, use.percentages = FALSE)
resultsDAC_NCvVU@results
View(resultsDAC_NCvVU@results)
write.csv(resultsDAC_NCvVU@results, "Naive_Chal_v_Vax_Unchal_DAC_p_values.csv", row.names = FALSE)

tiff("Naive_Chal vs Vax_Unchal.tiff",
     width=2000,
     height=1500, 
     res = 300)
SPADEVizR::plot(resultsDAC_NCvVU)
dev.off()


setwd(PrimaryDirectory)



###################################################
# Analysis with edgeR 
###################################################
edgeR_analysis <- function(data_experiment, data_control, experiment_sample_size, control_sample_size, export_file_name){
  
  # data_experiment = Vax_Chal.abundances
  #  data_control = Naive_Chal.abundances
  #  experiment_sample_size = 8
  #  control_sample_size = 7
  #  export_file_name = "Vax_Chal v Naive_Chal EdgeR Analysis"

  tmp_cluster.abundances = cbind(data_experiment, data_control)
  # Group assignment
  conditions1 = rep("A", experiment_sample_size)
  conditions2 = rep("B", control_sample_size)
  conditions = c(conditions1, conditions2)
  
  y <- DGEList(tmp_cluster.abundances)
  
  ## ------------------------------------------------------------------------
  #Because data has been rescaled, cant 
  #keep <- aveLogCPM(y) >= aveLogCPM(5, mean(y$samples$lib.size))
  #y <- y[keep,]
  
  y = y
  
  ## ------------------------------------------------------------------------
  design <- model.matrix(~factor(conditions))
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust=TRUE)
  res <- glmQLFTest(fit, coef=2)
  DAC = topTags(res, n=200, adjust.method="BH", sort.by="PValue", p.value=1)
  print(DAC)
  View(DAC)
  write.csv(DAC, paste(export_file_name, ".csv", sep = ""),row.names = TRUE)
}

Vax_Chal.abundances = cluster.abundances[, colnames(cluster.abundances) %in% Vax_Chal]
Naive_Chal.abundances = cluster.abundances[, colnames(cluster.abundances) %in% Naive_Chal]
Control.abundances = cluster.abundances[, colnames(cluster.abundances) %in% Control]
Vax_Unchal.abundances = cluster.abundances[, colnames(cluster.abundances) %in% Vax_Unchal]

dir.create("EdgeR Analysis", showWarnings = FALSE)
setwd("EdgeR Analysis")

### Change data_experiment and data_control to sample names you want to compare as well as 
### experiment_sample_size and control_sample_size to number of conditions in each sample
### finally, change export_file_name to depict groups being compared
edgeR_analysis(data_experiment = Vax_Chal.abundances,
               data_control = Naive_Chal.abundances,
               experiment_sample_size = 8,
               control_sample_size = 7,
               export_file_name = "Vax_Chal v Naive_Chal EdgeR Analysis")

edgeR_analysis(data_experiment = Vax_Chal.abundances,
               data_control = Vax_Unchal.abundances,
               experiment_sample_size = 8,
               control_sample_size = 8,
               export_file_name = "Vax_Chal v Vax_Unchal EdgeR Analysis")

edgeR_analysis(data_experiment = Control.abundances,
               data_control = Naive_Chal.abundances,
               experiment_sample_size = 8,
               control_sample_size = 7,
               export_file_name = "Control v Naive_Chal EdgeR Analysis")


edgeR_analysis(data_experiment = Control.abundances,
               data_control = Vax_Unchal.abundances,
               experiment_sample_size = 8,
               control_sample_size = 8,
               export_file_name = "Control v Vax_Unchal EdgeR Analysis")

setwd(PrimaryDirectory)

####################################
#SCATTER PLOT GENERATOR
###################################
Grouped_SheetFour <- read_excel("./Abundance DR3 Day 52 B Cell Data.xlsx", sheet = "Sheet 1")
write.table(Grouped_SheetFour, file = "Data for Scatter Plot Grouped.txt", sep = "\t",row.names = FALSE, col.names = TRUE)


#load data
data <- read.table("Data for Scatter Plot Grouped.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#data$Cluster <- gsub(" ", "_", data$Cluster, fixed = TRUE)
rownames(data) <- data$Cluster
data <- data[,-1]

#transpose the data for ploting
data <- t(data)
data <- as.data.frame(data)

#group assignment
group_data <- read.table("group assignment for grouped.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
group_data$sample <- trim(group_data$sample)
group_data$sample = gsub(" ", ".", group_data$sample, fixed = TRUE)
data$group <- group_data$group#[match(rownames(data), group_data$sample)]
data <- data[, c(ncol(data), 1:(ncol(data)-1))]

dir.create("Grouped_Scatterplots", showWarnings = FALSE)
setwd("Grouped_Scatterplots")
x_order = factor(data$group, levels=c("Naive Unchallenged", "Naive Challenged", "Vax Unchallenged","Vax Challenged"), ordered=TRUE)
for(i in 2:ncol(data)){
  scatter_plot <- 
    ggplot(data, aes_string(x = x_order, fill = "group", y = colnames(data)[i]))+ 
    geom_dotplot(binaxis = "y", stackdir = "centerwhole") +
    stat_summary(fun.y = "median", size=0.5, geom = 'line', aes(group=1))+
    stat_summary(
      fun.ymin = function(z) { quantile(z,0.25) },
      fun.ymax = function(z) { quantile(z,0.75) },
      fun.y = median,
      width = 0.2,
      geom = "errorbar") + 
    theme(axis.text.x = element_text(size = 25, face = "bold", vjust = 1.0, hjust = 1.0, angle = 45)) +
    theme(axis.text.y = element_text(size = 20, face = "bold", vjust = 0.5, hjust = 0.5, angle = 0)) +
    theme(legend.position = "none")
  ggsave(scatter_plot,
         width = 20,
         height = 15,
         dpi = 300,
         filename = paste(colnames(data)[i], ".jpeg", sep = ""))
}
setwd(PrimaryDirectory)



### Displays an heatmap representation summarizing phenotypes for the overall dataset
heatmapViewer(results)
####################################################################################################################################################################################
####################################################################################################################################################################################