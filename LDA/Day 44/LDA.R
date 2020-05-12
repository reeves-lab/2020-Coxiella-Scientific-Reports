#Function for Running LDA
#Joshua Hess

#Linear Discriminant Analysis function
RunLDA = function(input_to_lda,plot_colors,cv_type,num_cores,day){
  
  #Run parallel for speed
  registerDoParallel(cores = num_cores)
  
  #Set up 'Group' column of data to be a factor
  input_to_lda$Group = as.factor(input_to_lda$Group)
  levels(input_to_lda$Group) <- make.names(levels(factor(input_to_lda$Group)))
  input_to_lda = input_to_lda %>%
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
  
  #Run lda
  set.seed(433)
  LDA <- caret::train(Group ~ . ,
                      input_to_lda,
                      method = 'lda',
                      trControl = custom)

  #Extract LDA coefficients
  CoefficientsLDA <- as.data.frame(LDA$finalModel$scaling)
  rownames(CoefficientsLDA) = stringr::str_replace_all(rownames(CoefficientsLDA),
                                                      "[`]", "")
  colnames(CoefficientsLDA)[1] = "Coefficients of Linear Discriminants"
  write.xlsx(CoefficientsLDA, paste(day,levels(input_to_lda$Group)[1],
                                   levels(input_to_lda$Group)[2],
                                   "Coeff.xlsx",sep = "_"), 
             row.names = TRUE)
  
  #Set up Coefficients waterfall plot data
  CoefficientsLDA = CoefficientsLDA[order(-CoefficientsLDA$`Coefficients of Linear Discriminants`),
                                  , drop = FALSE]
  
  # create new column for names
  CoefficientsLDA$Group <- ifelse(CoefficientsLDA$`Coefficients of Linear Discriminants` < 0,
                                 levels(input_to_lda$Group)[1],
                                 levels(input_to_lda$Group)[2])
  
  #Plot coefficients in waterfall plot
  tmp_values = CoefficientsLDA$`Coefficients of Linear Discriminants`
  features = rownames(CoefficientsLDA)
  features = factor(features, levels = features[order(tmp_values)])
  axval = (max(abs(tmp_values)) + 0.03)
  
  plot_coef = ggplot(CoefficientsLDA, aes(x=features, y=tmp_values)) + 
    geom_bar(stat='identity', aes(fill=Group), width=.6)  +
    scale_fill_manual(name="Group", 
                      labels = c(levels(input_to_lda$Group)[1],levels(input_to_lda$Group)[2]), 
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
    labs(x = "Phenotypic Class", y = expression(paste("LDA Coeffiecient"))) + 
    scale_y_continuous(limits = c(-axval, axval)) +
    coord_flip()
  ggsave(paste(day,levels(input_to_lda$Group)[1],
               levels(input_to_lda$Group)[2],
               "coeff_plot.jpeg",sep="_"),plot = plot_coef, units = "mm", width = 50, height = 50)
  
}

############################
#Run ML models function to return top models based on accuracy
LDA_run = function (x){
  
  #Set up control parameters
  set.seed(123121)
  custom <- caret::trainControl(method = "boot",
                                number = 1000,
                                verboseIter = TRUE, 
                                classProbs = TRUE,
                                summaryFunction = multiClassSummary,
                                savePredictions="final",returnData = TRUE,returnResamp = "all")
  
  
  #Run the LDA model
  set.seed(123121)
  LDA_model <- caret::train(Group ~ . ,
                            x,
                            method = 'lda',
                            trControl = custom,
                            metric = "Accuracy")
  
  #Return the bayes model results
  return(LDA_model)
}


#Function for iterative classification
RunIterativeLDA = function(df,name_of_response,export_results,num_cores){
  
  #df: Dataframe containing the model predictors as columns and observations as rows
  #name_of_response: character. Name of the column in df that is your response variable
  #range_of_vars: Maximum number of predictors for your final model (only here for saving computational time)
  #export_results: Logical. If true, a results table of all models is exported.
  
  #Register parallel processing
  registerDoParallel(cores = num_cores)
  
  
  name_of_response = "Cat"
  df = input_LDA
  
  #Obtain the temporary set of names to permute (Make sure that the factor variable used to classify is last column)
  Group = df[,which(names(df)==name_of_response)]
  tmp_data = df[,-which(names(df)==name_of_response)]
  #Form a list containing the set of combinations
  comb_list = lapply(seq_len(4), FUN = function(x)combn(colnames(tmp_data), x))
  
  #Loop for running all of the models that you indicated
  result_list = list()
  for (k in 1:ncol(tmp_data)){
    
    k=1
    
    #Obtain only the combinations from the data in the list element
    tmp_subset = as.data.frame(comb_list[[k]])
    colnames(tmp_subset) = colnames(tmp_data)[k]
    tmp_subset = cbind(tmp_subset,Group)
    
    #Set up the response variable to be a factor with levels
    tmp_subset[,"Group"] = as.factor(tmp_subset[,"Group"])
    levels(tmp_subset[,"Group"]) <- make.names(levels(factor(tmp_subset[,"Group"])))
    
    #Run the classifiers over each column
    tmp_model = LDA_run(x = tmp_subset)
    
    #Add the list_per_range to the result_list so we can store all ranges in a single structure
    result_list[[k]] = tmp_model
  }
  
  #Get the names for your list object
  names(result_list) = colnames(tmp_data)
  
  return_list = list()
  for (j in 1:length(result_list)){
    
    #Extract hold out statistics from Confusion Matrix
    tmp_confusion = confusionMatrix(result_list[[j]]$pred$pred,
                                    result_list[[j]]$pred$obs,mode = "everything",
                                    positive = result_list[[j]]$levels[2])
    tmp_confusion_frame = cbind(as.data.frame(t(tmp_confusion$overall)),as.data.frame(t(tmp_confusion$byClass)))
    colnames(tmp_confusion_frame) = paste("Resamples_",colnames(tmp_confusion_frame),sep = "")
    
    #Add the best model to the bayes_export_table dataframe
    tmp_confusion_frame$Index = paste("result_list","[[",j,"]]","$bestTune",sep = "")
    tmp_confusion_frame$Model = paste(as.character(result_list[[j]]$coefnames),collapse = "")
    tmp_confusion_frame$Model = gsub("`", "", tmp_confusion_frame$Model)
    
    #Ensure that bestTune is the model we chose!!!
    tmp_confusion_frame$bestTune_validation = paste(result_list[[j]]$bestTune,collapse = " , ")
    
    #Get the coefficient from the LDA model
    tmp_confusion_frame$Coefficient = result_list[[j]]$finalModel$scaling
    
    #Combine resampled stats with the full stats model
    final_results_frame = as.data.frame(tmp_confusion_frame)
    
    #add the model to the list
    return_list[[j]] = final_results_frame
  }
  return_table = do.call(rbind,return_list)
  if (export_results){
    write.xlsx(return_table,'Individual_LDA_results.xlsx')
  }
  #Return the table 
  return(return_table)
}






