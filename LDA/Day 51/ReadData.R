#Script for conducting LDA analysis in R
#Joshua Hess


#Check for missing packages and install if needed
list.of.packages <- c("caret","dplyr","openxlsx","magrittr","doParallel","MLmetrics")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
library(dplyr)
library(openxlsx)
library(caret)
library(doParallel)
library(MLmetrics)



ReadData = function(filename,sheetname,scale=TRUE,cluster_prefix=NULL){
  #This function will do the calculation for your correlation analysis. This analysis is made for data that
  #has cluster data attached to columns that correspond to the measures that you want to correlate to. The column
  #names in your excel file must be written as "Cluster #" for each cluster. See example excel file if you are having
  #format issues.
  
  
  #filename: Indicates the path to your excel file (Ex: filename = "user/path/../myexcelfile.xlsx"). You can input a list of paths as well in order 
  #to concatenate a list of files for correlation analysis (Ex: filename = list()). If you choose to do this, you must also include a sheet list 
  
  #export: Logical indicating whether or not you want to export excel files for the correlation data (includes p-value and R value)
  
  #cluster_prefix: List indicating what prefixes you want to add to each of your cluster types in the dataframe. You only need to input this if you are
  #including multiple filenames (Ex: cluster_prefix = list("T","B","I"))
  
  #Read the excel file(s):
  if (class(filename) == "list"){
    print("Detected a list of data files...")
    #Get the last day group assignment
    group = na.omit(read.xlsx(filename[[length(filename)]],sheet=sheetname[[length(sheetname)]]))
    group = dplyr::select(group, contains("Category"))
    colnames(group)="Cat"
    fin_group = subset(group, Cat == "Naive Challenge" | Cat == "Vaccinated Challenge")
    mice=na.omit(read.xlsx(filename[[length(filename)]],sheet=sheetname[[length(sheetname)]]))
    mice = cbind(mice,group)
    mice = subset(mice, Cat == "Naive Challenge" | Cat == "Vaccinated Challenge") %>%
      dplyr::select(contains("X"))
    #Get the data
    data = list()
    clusters = list()
    not_clusters = list()
    for (i in 1:length(filename)){
      data[[i]] = na.omit(read.xlsx(filename[[i]],sheet=sheetname[[i]]))
      #Scale the data if chosen
      clusters[[i]] = dplyr::select(data[[i]], starts_with("Cluster"))
      if (scale){
        clusters[[i]] = scale(clusters[[i]],center = TRUE,scale = TRUE)
      }
      data[[i]] = cbind(data[[i]],group)
      data[[i]] = subset(data[[i]], Cat == "Naive Challenge" | Cat == "Vaccinated Challenge")
      clusters[[i]] = cbind(clusters[[i]],group)
      clusters[[i]] = subset(clusters[[i]], Cat == "Naive Challenge" | Cat == "Vaccinated Challenge") %>%
        dplyr::select(-contains("Cat"))
      colnames(clusters[[i]]) = paste(cluster_prefix[[i]],colnames(clusters[[i]]),sep = " ")
      not_clusters[[i]] = dplyr::select(data[[i]],-contains("Cluster")) %>%
        dplyr::select(-contains("Category")) %>%
        dplyr::select(-contains("Cat")) %>%
        dplyr::select(-contains("X"))
    }
    clusters = do.call(cbind,clusters)
    not_clusters = do.call(cbind,not_clusters)
    #Remove duplicate columns
    not_clusters = not_clusters[,!duplicated(colnames(not_clusters),fromLast = TRUE)] %>%
      dplyr::select(-contains("Total"))
  }
  #Create a list object from the cluster data and the clinical data
  data = list(clusters,not_clusters,mice,fin_group);names(data) = c("Clusters","Clinical","ID","Group")
  #Return the data
  return(data)
}


