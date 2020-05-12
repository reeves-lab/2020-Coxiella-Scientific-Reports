
#Function for iterative Naive Bayes Classification
run.ML.iterations = function(df,name_of_response,range_of_vars,export_results,num_cores){
  
  #df: Dataframe containing the model predictors as columns and observations as rows
  #name_of_response: character. Name of the column in df that is your response variable
  #range_of_vars: Maximum number of predictors for your final model (only here for saving computational time)
  #export_results: Logical. If true, a results table of all models is exported.
  
  #Register parallel processing
  registerDoParallel(cores = num_cores)
  
  #Run ML models function to return top models based on accuracy
  nested.ML.run = function (x){
    
    #Extract the combinations from a column in the data frame tmp_subset
    tmp_names = append(as.character(x),name_of_response, after = length(as.character(x)))
    tmp_input = df[,which(names(df) %in% tmp_names)]
    
    #Set up the response variable to be a factor with levels
    tmp_input[,ncol(tmp_input)] = as.factor(tmp_input[,ncol(tmp_input)])
    levels(tmp_input[,ncol(tmp_input)]) <- make.names(levels(factor(tmp_input[,ncol(tmp_input)])))
    
    #Set up control parameters
    set.seed(123121)
    custom <- caret::trainControl(method = "boot",
                                  number = 1000,
                                  verboseIter = TRUE, 
                                  classProbs = TRUE,
                                  summaryFunction = multiClassSummary,
                                  savePredictions="final",returnData = TRUE,returnResamp = "all")
    
    #Run the naive bayes model
    bayes_grid <- expand.grid(
      usekernel = c(TRUE, FALSE),
      adjust = seq(0, 2, by = 1),
      fL = 0:2)
    
    set.seed(123121)
    Bayes_model <- caret::train(Group ~ . ,
                                tmp_input,
                                method = 'nb',
                                trControl = custom,
                                tuneGrid = bayes_grid,
                                metric = "Accuracy")
    
    #Run the LDA model
    set.seed(123121)
    LDA_model <- caret::train(Group ~ . ,
                              tmp_input,
                              method = 'lda',
                              trControl = custom,
                              metric = "Accuracy")
    
    #Run the Linear SVM model
    set.seed(123121)
    svm_grid = expand.grid(C = seq(1,10,by = 0.5))
    SVM_model <- caret::train(Group ~., 
                              tmp_input, 
                              method = "svmLinear",
                              trControl = custom,
                              tuneGrid = svm_grid,
                              metric = "Accuracy")
    
    #Run the neural network model
    nnet_grid <- expand.grid(decay = c(0.5,0.1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7), 
                             size = c(3,5,10,12,15,20,25))
    
    set.seed(123121)
    nnet_model <- caret::train(Group ~., 
                               tmp_input, 
                               method = "nnet",
                               trControl = custom,
                               tuneGrid = nnet_grid,
                               metric = "Accuracy")
    
    #Return the bayes model results
    return(list(Bayes_model,LDA_model,SVM_model,nnet_model))
  }
  
  #Obtain the temporary set of names to permute (Make sure that the factor variable used to classify is last column)
  tmp_data = df[,-which(names(df)==name_of_response)]
  
  #Form a list containing the set of combinations
  comb_list = lapply(seq_len(ncol(tmp_data)), FUN = function(x)combn(colnames(tmp_data), x)) 
  
  #Loop for running all of the models that you indicated
  result_list = list()
  for (k in range_of_vars){
    
    #Obtain only the combinations from the data in the list element
    tmp_subset = as.data.frame(comb_list[[k]]) 
    
    #Run the classifiers over each column
    tmp_list = apply(tmp_subset, 2, nested.ML.run)
    
    #Add the list_per_range to the result_list so we can store all ranges in a single structure
    result_list[[k]] = tmp_list
  }
  
  
  #Extract the highest accuracy model by cross validation from each result in result_list
  bayes_export_list = list(0)
  LDA_export_list = list(0)
  SVM_export_list = list(0)
  nnet_export_list = list(0)
  t=1
  for (j in 1:length(result_list)){
    for (k in 1:length(result_list[[j]])){
      for (h in 1:length(result_list[[j]][[k]])){
        
        #Extract hold out statistics from Confusion Matrix
        tmp_confusion = confusionMatrix(result_list[[j]][[k]][[h]]$pred$pred,
                                        result_list[[j]][[k]][[h]]$pred$obs,mode = "everything",
                                        positive = result_list[[j]][[k]][[h]]$levels[2])
        tmp_confusion_frame = cbind(as.data.frame(t(tmp_confusion$overall)),as.data.frame(t(tmp_confusion$byClass)))
        colnames(tmp_confusion_frame) = paste("Resamples_",colnames(tmp_confusion_frame),sep = "")
        
        #Add the best model to the bayes_export_table dataframe
        tmp_confusion_frame$Index = paste("result_list","[[",j,"]]","[[",k,"]]","[[",h,"]]","$bestTune",sep = "")
        tmp_confusion_frame$Model = paste(as.character(result_list[[j]][[k]][[h]]$coefnames),collapse = "")
        
        #Ensure that bestTune is the model we chose!!!
        tmp_confusion_frame$bestTune_validation = paste(result_list[[j]][[k]][[h]]$bestTune,collapse = " , ")
        
        #Combine resampled stats with the full stats model
        final_results_frame = as.data.frame(tmp_confusion_frame)
        
        #Take only the first model to keep resuts clean as a table (will be the same model chosen by caret 'bestTune')
        if (h==1){
          bayes_export_list[[t]] = final_results_frame[1,]
        }else if (h==2){
          LDA_export_list[[t]] = final_results_frame[1,]
        }else if (h==3){
          SVM_export_list[[t]] = final_results_frame[1,]
        }else if (h==4){
          nnet_export_list[[t]] = final_results_frame[1,]
        }
        t = t+1
      }
    }
  }
  bayes_export_table = do.call(rbind,bayes_export_list)
  LDA_export_table = do.call(rbind,LDA_export_list)
  SVM_export_table = do.call(rbind,SVM_export_list)
  nnet_export_table = do.call(rbind,nnet_export_list)
  
  #Logical to export the results talbe
  if (export_results ==TRUE){
    write.xlsx(bayes_export_table,paste("Naive_Bayes_results.xlsx"))
    write.xlsx(LDA_export_table,paste("LDA_results.xlsx"))
    write.xlsx(SVM_export_table,paste("SVM_results.xlsx"))
    write.xlsx(nnet_export_table,paste("Nnet_results.xlsx"))
  }
  
  #Return the result_list so that you can extract any model you wish
  return(list(result_list,bayes_export_table,LDA_export_table,SVM_export_table,nnet_export_table))
}




