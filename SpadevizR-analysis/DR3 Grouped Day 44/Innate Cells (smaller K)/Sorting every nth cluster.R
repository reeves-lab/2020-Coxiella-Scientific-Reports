
data <- read.table("DR3 D51 T Cells Pearsons Coefficient.txt", header = T, sep = "\t", stringsAsFactors = F)

data <- data[order(data$experiment1_cluster, data$cor, decreasing = TRUE),]

cluster.name <- unique(data$experiment1_cluster)

new.data <- data
new.data <- new.data[1,]
new.data[1, ] <- NA


for(i in cluster.name){
  #i = "Grp1_Cluster_5169"
  tmp.data <- data[which(data$experiment1_cluster == i), ]
  tmp.data <- tmp.data[1:5, ]
  new.data <- rbind(new.data, tmp.data)
}

new.data <- na.omit(new.data)

write.table(new.data, "DR3 D51 T Cells smaller K Sorted Pearsons Data.txt", sep = "\t", row.names = F)
