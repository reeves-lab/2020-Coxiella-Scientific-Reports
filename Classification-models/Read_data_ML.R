#Function for reading ML models data for each day
read.data = function(excel_file,cell_type,excel_sheet){
  
  #Read the excel files
  tmp_data = read.xlsx(paste(excel_file),sheet = excel_sheet)
  
  #na.omit just to be sure
  tmp_data <- na.omit(tmp_data)
  
  #Extract the Category and assign globally
  Group <<- tmp_data[c("Category")]
  
  #Get only the cluster data from the datasheet
  cluster_data <- dplyr::select(tmp_data, starts_with("Cluster"))
  
  #Set the rownames to be nice
  rownames(cluster_data)[1:16] <- paste("Grp1 Mouse", seq(1,16), sep = " ")
  rownames(cluster_data)[17:31] <- paste("Grp2 Mouse", seq(2,16), sep = " ")
  
  #Set the column names to be nice
  if (excel_sheet == "Day 10"){
    colnames(cluster_data)[1:ncol(cluster_data)] = paste(cell_type, "D10 Cluster", 
                                                         suffix = seq(1:ncol(cluster_data)), sep = " ")
  }else if (excel_sheet == "Day 35"){
    colnames(cluster_data)[1:ncol(cluster_data)] = paste(cell_type, "D35 Cluster", 
                                                         suffix = seq(1:ncol(cluster_data)), sep = " ")
  }else if (excel_sheet == "Day 44"){
    colnames(cluster_data)[1:ncol(cluster_data)] = paste(cell_type, "D44 Cluster", 
                                                         suffix = seq(1:ncol(cluster_data)), sep = " ")
  }else if (excel_sheet == "Day 51"){
    colnames(cluster_data)[1:ncol(cluster_data)] = paste(cell_type, "D51 Cluster", 
                                                         suffix = seq(1:ncol(cluster_data)), sep = " ")
  }else if (excel_sheet == "Day 52"){
    colnames(cluster_data)[1:ncol(cluster_data)] = paste(cell_type, "D52 Cluster", 
                                                         suffix = seq(1:ncol(cluster_data)), sep = " ")
  }
  
  #Return scaled cluster data so that we can use it to cbind
  return(as.data.frame(scale(cluster_data,center = TRUE,scale = TRUE)))
}