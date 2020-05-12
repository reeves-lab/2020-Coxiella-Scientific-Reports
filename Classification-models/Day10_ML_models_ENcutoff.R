#Day 10 Classification iterations
#Joshua Hess

library(openxlsx)
library(caret)
library(ggplot2)
library(dplyr)
library(stringr)
library(magrittr)
library(varhandle)
library(doParallel)
library(e1071)
library(klaR)
library(MLmetrics)
library(nnet)
source("Read_data_ML.R")
source("ML_models_exhaustive.R")

#Read Day 10 Data
cluster_T = read.data(excel_file = "HistoScore and Phenotype Data_T Cell Clusters(smaller K)_DR3.xlsx",
                      cell_type = "T",
                      excel_sheet = "Day 10")
cluster_B = read.data(excel_file = "HistoScore and Phenotype Data_B Cell Clusters(smallerK)_DR3.xlsx",
                      cell_type = "B",
                      excel_sheet = "Day 10")
cluster_I = read.data(excel_file = "HistoScore and Phenotype Data_Innate Cell Clusters(mixed K)_DR3.xlsx",
                      cell_type = "I",
                      excel_sheet = "Day 10")

data_full = cbind(cluster_T,
                  cluster_B,
                  cluster_I)

#Change colnames in Group dataframe so that we can input to the classifiers
colnames(Group) = "Group"

#Read Prism data and include only those clusters that show up in top 10 elastic net
data_EN = read.xlsx("D10_Prism_EN_Coefficients_Log Odds.xlsx",sheet = "Prism Data")
data_EN <-data_EN[abs(data_EN$Change.in.Log.Odds.per.Unit.Increase) >= 0.02,] %>% 
  top_n(10, wt = abs(Change.in.Log.Odds.per.Unit.Increase)) %>%
  arrange(desc(abs(Change.in.Log.Odds.per.Unit.Increase))) %>%
  as.data.frame()

tmp_input = data_full[,which(colnames(data_full) %in% c(data_EN$X1))]
data_input = cbind(tmp_input,Group)


#Run the ML classifiers
ML_results_list = run.ML.iterations(df = data_input,
                                    name_of_response = c("Group"),
                                    range_of_vars = 1:10,
                                    export_results = TRUE,
                                    num_cores = 28)


















