#Do all spearman vs pearson comparisons
rescale_to_0_1 <- function(experiment_name, experiment_file, rescale = TRUE){
  #This function will resscale your data from 0 to 1
  
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

Run.Correlations = function(experiment1_file,experiment2_file,export_file){
  #Function for exporting the correlation results
  
  #Read the data and rescale
  experiment1 = rescale_to_0_1(experiment_name = "Grp1",
                               experiment_file = experiment1_file,
                               rescale = TRUE)
  
  experiment2 = rescale_to_0_1(experiment_name = "Grp2",
                               experiment_file = experiment2_file,
                               rescale = TRUE)

  
  experiment1_1 = experiment1
  experiment1 = experiment1[, colnames(experiment1) != "Count1"]
  
  experiment2_1 = experiment2
  experiment2 = experiment2[, colnames(experiment2) != "Count1"]
  
  #create a blank table to store the pearson correlation results
  corr_data<-data.frame(Experiment1_Cluster = 0, Experiment2_Cluster = 0, Experiment1_Count = 0, Experiment2_Count = 0)
  
  #perform pairwise correlations between experiment1 and experiment2
  t=1
  for(i in 1:nrow(experiment1)){
    for(j in 1:nrow(experiment2)){
      corr_data$Experiment1_Cluster[t]<-rownames(experiment1)[i]
      corr_data$Experiment2_Cluster[t]<-rownames(experiment2)[j]
      
      pearson_stats<-cor.test(as.numeric(experiment1[i,]),as.numeric(experiment2[j,]),method = "pearson")
      spearman_stats<-cor.test(as.numeric(experiment1[i,]),as.numeric(experiment2[j,]),method = "spearman")
      
      corr_data$`Pearsons Coefficient`[t]<-pearson_stats$estimate
      corr_data$`Spearmans Coefficient`[t]<-spearman_stats$estimate
      corr_data$`Pearsons p-value`[t]<-pearson_stats$p.value
      corr_data$`Spearmans p-value`[t]<-spearman_stats$p.value
      
      corr_data$Experiment1_Count[t]<-experiment1_1$Count1[i]
      corr_data$Experiment2_Count[t]<-experiment2_1$Count1[j]
      
      t<-t+1
      corr_data<-rbind(corr_data, 0)
    }
  }
  
  #Do lots of data preprocessing with long names
  data_match = read.table("Matching Clusters for Visualization.txt",sep = "\t");data_match = data_match[,colnames(data_match) %in% c("V1","V2","V3")]
  data_match <- data.frame(lapply(data_match, as.character), stringsAsFactors=FALSE);colnames(data_match)=data_match[1,]
  data_match = data_match[-1,];data_match$Experiment1_Cluster = paste("Grp1_Cluster_",data_match$Experiment1_Cluster,sep = "")
  data_match$Experiment2_Cluster=paste("Grp2_Cluster_",data_match$Experiment2_Cluster,sep = ""); colnames(data_match)[3] = "Matched Cluster #"
  
  #Only extract the matched clusters from the correlation data
  final = merge(corr_data, data_match,by=c("Experiment1_Cluster","Experiment2_Cluster"))
  #Sorting the data for better view
  final = final[order(as.numeric(final$`Matched Cluster #`), decreasing = FALSE),]
  
  # create a CSV file storing the pearson correlation data
  write.csv(final, export_file, row.names = FALSE)  
}
