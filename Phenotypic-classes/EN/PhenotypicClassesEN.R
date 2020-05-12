

#Phenotypic Families Elastic Net Regression
library(openxlsx)
library(caret)
library(ggplot2)
library(dplyr)
library(stringr)
library(magrittr)
library(varhandle)
library(doParallel)

#Function for reading Day 10 and Day 35 data
read.data.D10.D35 = function(cell_list,days){
  
  # days = "D10"
  # cell_list = Cell_types[2]

  #read data for each cell type in this function
  temp_data = read.xlsx("Elastic Net input_phenotypic_families_Updated Bcells_20phenoFams.xlsx",
                        sheet = paste(cell_list, "cells D10_D35",
                                      sep = " "))
  
  #Append data to be in format for taking median of each family abundance
  colnames(temp_data) = paste(colnames(temp_data),temp_data[1,],sep = " ")
  temp_data = temp_data[grep(days, temp_data$Category), ]
  
  temp_data = temp_data[,!grepl("c01_Grp2",names(temp_data))] %>%
    select(-c("Category Cluster type day number")) %>%
    t()%>%
    as.data.frame()
    
  #Z score each column so that relative abundances are on the same scale for each family
  colnames(temp_data) <- paste(cell_list,as.character(unlist(temp_data[1,])),sep = " ")
  temp_data_pip = temp_data[!rownames(temp_data) %in% c("Family NA"), ]
  row_names_scaled = rownames(temp_data_pip)
  
  #Convert data to numeric so that we can scale
  temp_data_pip = unfactor(temp_data_pip)
  temp_data_pip = as.data.frame(sapply(temp_data_pip, as.numeric))
  
  #Z score data 
  temp_data_scaled = scale(temp_data_pip, scale = TRUE,center = TRUE)
  
  #Create column for family assignment
  family_assign = colnames(temp_data_scaled)
  rownames(temp_data_scaled) = row_names_scaled
  families_table = as.data.frame(t(temp_data_scaled))
  families_table$Family = family_assign
  
  #Create median stats table for frequency of each family in mouse
  stats_table = do.call(data.frame,
                        aggregate(families_table[,1:ncol(families_table)-1],
                        by = list(families_table$Family),
                        function(x) c(median = median(x))))
  
}

#Read frequency data for phenotypic families
Cell_types = c("T","B","I")
day_10_list = lapply(Cell_types,function(x){read.data.D10.D35(x,c("D10"))})
day_35_list = lapply(Cell_types,function(x){read.data.D10.D35(x,c("D35"))})


#Set up data for the elastic net regression
data.cleanup.D10.D35 = function(x){
  
  #Combine all data from the list
  data = t(do.call(rbind,x))
  colnames(data) = data[1,]
  data_fin = as.data.frame(data[-1,])
  
  #Set group assignment for each row (Naive vs Vaccinated)
  data_fin$Group[grepl("Naive",rownames(data_fin))] = "Naive"
  data_fin$Group[grepl("Vaccinated",rownames(data_fin))] = "Vaccinated"
  return(as.data.frame(data_fin))
}

#Clean up day 10 and day 35 data
day10_input = data.cleanup.D10.D35(day_10_list)
day35_input = data.cleanup.D10.D35(day_35_list)


#Function for reading Day 44 and Day 51 data
read.data.D44.D51 = function(cell_list,days){
  
  #read data for each cell type in this function
  temp_data = read.xlsx("Elastic Net input_phenotypic_families_Updated Bcells_20phenoFams.xlsx",
                        sheet = paste(cell_list, "cells D44_D51",
                                      sep = " "))
  
  #Append data to be in format for taking median of each family abundance
  # colnames(temp_data) = paste(colnames(temp_data),temp_data[1,],sep = " ")
  temp_data = temp_data[grep(days, temp_data$Category), ]
  colnames(temp_data)[1] = "Family NA"
  temp_data = t(temp_data) %>%
    as.data.frame()
  
  #Z score each column so that relative abundances are on the same scale for each family
  colnames(temp_data) <- paste(cell_list,as.character(unlist(temp_data[1,])),sep = " ")
  temp_data_pip = temp_data[!rownames(temp_data) %in% c("Family.NA","Category"), ]
  row_names_scaled = rownames(temp_data_pip)
  
  #Convert data to numeric so that we can scale
  temp_data_pip = unfactor(temp_data_pip)
  temp_data_pip = as.data.frame(sapply(temp_data_pip, as.numeric))
  
  #Z score data 
  temp_data_scaled = scale(temp_data_pip, scale = TRUE,center = TRUE)
  
  #Create column for family assignment
  family_assign = colnames(temp_data_scaled)
  rownames(temp_data_scaled) = row_names_scaled
  families_table = as.data.frame(t(temp_data_scaled))
  families_table$Family = family_assign
  
  #Create median stats table for frequency of each family in mouse
  stats_table = do.call(data.frame,
                        aggregate(families_table[,1:ncol(families_table)-1],
                                  by = list(families_table$Family),
                                  function(x) c(median = median(x))))
  
}

day_44_list = lapply(Cell_types,function(x){read.data.D44.D51(x,c("D44"))})
day_51_list = lapply(Cell_types,function(x){read.data.D44.D51(x,c("D51"))})

#Set up data for the elastic net regression
data.cleanup.D44.D51 = function(x){
  
  #Combine all data from the list
  data = t(do.call(rbind,x))
  colnames(data) = data[1,]
  data_fin = as.data.frame(data[-1,])
  
  #Set group assignment for each row (Naive vs Vaccinated)
  data_fin$Group[grepl("Naive.Challenge",rownames(data_fin))] = "Naive Challenge"
  data_fin$Group[grepl("Vaccinated.Challenge",rownames(data_fin))] = "Vaccinated Challenge"
  return(as.data.frame(data_fin))
}

day44_input = data.cleanup.D44.D51(day_44_list) %>%
  na.omit()
day51_input = data.cleanup.D44.D51(day_51_list) %>%
  na.omit()

#---------- Elastic Net regression -------------------

#Elastic Net function
run.elastic.net = function(input_to_reg,plot_colors,cv_type,num_cores,day){

  #Run parallel for speed
  registerDoParallel(cores = num_cores)
  
  #z score predictor variables so that regualrization is equally weighted....data already scaled here
  # input_to_reg = input_to_reg %>%
  #   mutate_at(vars(-Group),funs(c(scale(.))))
  
  #Set up 'Group' column of regression data to be a factor
  input_to_reg$Group = as.factor(input_to_reg$Group)
  levels(input_to_reg$Group) <- make.names(levels(factor(input_to_reg$Group)))
  input_to_reg = input_to_reg %>%
    mutate_at(vars(-Group),funs(c(unfactor(.))))
  
  #Set up cross validation parameters - bootstrap CV or repeated CV
  if (cv_type == "bootstrap"){
    custom <- caret::trainControl(method = "boot",
                                  number = 1000,
                                  p = 0.75,
                                  verboseIter = TRUE, 
                                  classProbs = TRUE, summaryFunction = twoClassSummary)
  } else if (cv_type == "repeated"){
    custom <- caret::trainControl(method = "repeatedcv",
                                  number = 4,
                                  repeats = 100,
                                  verboseIter = TRUE, classProbs = TRUE)
  }
  
  #Run logistic elastic net regression
  set.seed(433)
  ElasticNet <- caret::train(Group ~ . ,
                             input_to_reg,
                             method = 'glmnet',family = "binomial", standardize = FALSE,
                             tuneGrid = expand.grid(
                               lambda = seq(0.00000001, 1,length = 100),
                               alpha = seq(0,1,length = 50)),
                             trControl = custom)
  #Plot cross validation
  ggplot(data = ElasticNet) +
    theme(text = element_text(size=17),         
          axis.text.x = element_text(face="bold", color="black", size=23, angle=0),   
          axis.text.y = element_text(face="bold", color="black", size=23, angle=0),    
          panel.background = element_blank(),                                         
          axis.line = element_line(colour = "black"),                                 
          panel.grid = element_line(colour = "grey"),                                 
          plot.margin = unit(c(2,2,2,2), "cm"),                                   
          axis.title.y = element_text(angle = 90, vjust = 10, hjust = 0.5, size = 25),
          axis.title.x = element_text(angle = 0, vjust = -10, hjust = 0.5, size = 25),
          legend.title = element_text(size=19), 
          legend.text = element_text(size=19))
  ggsave(paste(day,levels(input_to_reg$Group)[1],
               levels(input_to_reg$Group)[2],"Cross_valid.jpeg",sep = "_"),
         plot = last_plot(), units = "mm",
         width = 500, height = 300, limitsize = FALSE)
  
  #Export cross validation summary
  cvSummary = as.data.frame(ElasticNet$results)
  write.xlsx(cvSummary, paste(day,levels(input_to_reg$Group)[1],
                              levels(input_to_reg$Group)[2],
                              "Cross_valid_summary.xlsx",sep = "_"))
  
  #Extract log odds coefficients
  bestEN <- ElasticNet$finalModel
  CoefficientsEN = as.matrix(coef(bestEN, s = ElasticNet$bestTune$lambda))
  CoefficientsEN = as.data.frame(CoefficientsEN)
  rownames(CoefficientsEN) = stringr::str_replace_all(rownames(CoefficientsEN),
                                                      "[`]", "")
  colnames(CoefficientsEN)[1] = "Change in Log Odds per Unit Increase"
  write.xlsx(CoefficientsEN, paste(day,levels(input_to_reg$Group)[1],
                                   levels(input_to_reg$Group)[2],
                                   "Coeff.xlsx",sep = "_"), 
             row.names = TRUE)
  
  #Set up Coefficients waterfall plot data
  CoefficientsEN = subset(CoefficientsEN, CoefficientsEN$`Change in Log Odds per Unit Increase`!=0)
  CoefficientsEN = CoefficientsEN[2:nrow(CoefficientsEN),,drop = FALSE]
  CoefficientsEN = CoefficientsEN[order(-CoefficientsEN$`Change in Log Odds per Unit Increase`),
                                  , drop = FALSE]
  
  CoefficientsEN$`Group` <- rownames(CoefficientsEN)  # create new column for car names
  CoefficientsEN$Group <- ifelse(CoefficientsEN$`Change in Log Odds per Unit Increase` < 0,
                                 levels(input_to_reg$Group)[1],
                                 levels(input_to_reg$Group)[2])
  
  #Plot coefficients in waterfall plot
  logoddsvalue = CoefficientsEN$`Change in Log Odds per Unit Increase`
  features = rownames(CoefficientsEN)
  features = factor(features, levels = features[order(logoddsvalue)])
  axval = (max(abs(logoddsvalue)) + 0.03)
  
  plot_coef = ggplot(CoefficientsEN, aes(x=features, y=logoddsvalue)) + 
    geom_bar(stat='identity', aes(fill=Group), width=.6)  +
    scale_fill_manual(name="Group", 
                      labels = c(levels(input_to_reg$Group)[1],levels(input_to_reg$Group)[2]), 
                      values = plot_colors
    ) + 
    labs(subtitle="", title= "") +
    theme(text = element_text(size=0.5),                         
          axis.text.x = element_text(face="bold", color="black", size=2.5, angle=0),
          axis.text.y = element_text(face="bold", color="black", size=2.5, angle=0),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size = 0.15),
          panel.grid = element_line(colour = "grey", size = 0.1),
          plot.margin = unit(c(2,2,2,2), "mm"),
          axis.title.y = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 4),
          axis.title.x = element_text(angle = 0, vjust = -0.3, hjust = 0.5, size = 4),
          plot.title = element_text(vjust=3, hjust = 0),
          plot.subtitle = element_text(vjust = 3),
          legend.title = element_text(size=2), 
          legend.text = element_text(size=2),
          legend.key.size = unit(0.3,"line"),
          legend.margin = ggplot2::margin(t=-0.5, unit = "mm")) +
    labs(x = "Phenotypic Class", y = expression(paste(Delta, "Log Odds/","Unit Increase"))) + 
    scale_y_continuous(limits = c(-axval, axval)) +
    coord_flip()
  ggsave(paste(day,levels(input_to_reg$Group)[1],
               levels(input_to_reg$Group)[2],
               "coeff_plot.jpeg",sep="_"),plot = plot_coef, units = "mm", width = 50, height = 50)
  
}
  
run.elastic.net(input_to_reg = day10_input,
                 plot_colors = c("black","grey60"),cv_type = "repeated",
                 num_cores = 15, day = "D10")
run.elastic.net(input_to_reg = day35_input,
                plot_colors = c("black","grey60"),cv_type = "repeated",
                num_cores = 15, day = "D35")
run.elastic.net(input_to_reg = day44_input,
                plot_colors = c("black","grey60"),cv_type = "repeated",
                num_cores = 15, day = "D44")
run.elastic.net(input_to_reg = day51_input,
                plot_colors = c("black","grey60"),cv_type = "repeated",
                num_cores = 15, day = "D51")


