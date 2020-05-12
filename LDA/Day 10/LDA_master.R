#Master script for LDA Analysis
#Joshua Hess

#Get the functions that you will need
source("ReadData.R")
source("LDA_models_comb.R")

#Read the data
filename = list("HistoScore and Phenotype Data_B Cell Clusters(smallerK)_DR3.xlsx",
                "HistoScore and Phenotype Data_B Cell Clusters(smallerK)_DR3.xlsx",
                "HistoScore and Phenotype Data_B Cell Clusters(smallerK)_DR3.xlsx",
                "HistoScore and Phenotype Data_B Cell Clusters(smallerK)_DR3.xlsx",
                "HistoScore and Phenotype Data_Innate Cell Clusters(mixed K)_DR3.xlsx",
                "HistoScore and Phenotype Data_Innate Cell Clusters(mixed K)_DR3.xlsx",
                "HistoScore and Phenotype Data_Innate Cell Clusters(mixed K)_DR3.xlsx",
                "HistoScore and Phenotype Data_Innate Cell Clusters(mixed K)_DR3.xlsx",
                "HistoScore and Phenotype Data_T Cell Clusters(smaller K)_DR3.xlsx",
                "HistoScore and Phenotype Data_T Cell Clusters(smaller K)_DR3.xlsx",
                "HistoScore and Phenotype Data_T Cell Clusters(smaller K)_DR3.xlsx",
                "HistoScore and Phenotype Data_T Cell Clusters(smaller K)_DR3.xlsx")

sheetname = list("Day 10","Day 35","Day 44","Day 51",
                 "Day 10","Day 35","Day 44","Day 51",
                 "Day 10","Day 35","Day 44","Day 51")

cluster_prefix = list("B D10","B D35","B D44","B D51",
                      "I D10","I D35","I D44","I D51",
                      "T D10","T D35","T D44","T D51")

data = ReadData(filename,sheetname,scale = TRUE,cluster_prefix)

#Get the data in the format the we need
tmp_frame = cbind(data[["Clusters"]],data[["ID"]],data[["Group"]])
clin = cbind(data[["ID"]],data[["Group"]],data[["Clinical"]])
#Get only the vaccinated mice that have no splenic burden
vax_sick = clin[clin$`Spleen.Burdon.(copies/mg)`!=0&clin$Cat=="Vaccinated Challenge",]$X
#Remove the sick mice from the tmp_frame
input_LDA = tmp_frame[-which(tmp_frame$X %in% vax_sick),] %>%
  dplyr::select(-contains("X"))%>%
  rename(Group = Cat)

#Read the filtered data
filt = read.xlsx("Individual_LDA_results Filtered for Day.xlsx",sheet = "Day10")%>%
  arrange(desc(Resamples_Accuracy))%>%
  dplyr::select(contains("Model"))%>%
  #Take top n rows ordered by Model
  slice(1:15)

filt_names = c(filt$Model,"Group")
filt_data = input_LDA[,which(names(input_LDA) %in% filt_names)]

#Run the calculations
Result_List = run.ML.iterations(df = filt_data,
                                name_of_response = c("Group"),
                                export_results = TRUE,
                                max_size = 3,
                                num_cores = 25)


#Extract the predictions and plots for all models
table_stats = list()
index = 1
setwd("LDA Projections")
for (i in 1:length(Result_List[[1]])){
  for (j in 1:length(Result_List[[1]][[i]])){
    for (k in 1:length(Result_List[[1]][[i]][[j]])){
      
      #Extract the model
      mod = Result_List[[1]][[i]][[j]][[k]]$finalModel
      clusters = gsub("`","",mod[["xNames"]])
      cols = c(clusters,"Cat")
      
      #Change the columns to be vax sick and vax healthy
      tmp_frame[which(tmp_frame$X %in% vax_sick),]$Cat = "Vaccinated Challenge Sick"
      
      #Insert the predictors that were used for this model for all test subjects
      mod_dat = tmp_frame[,which(names(tmp_frame)%in%cols)]
      
      #Get the cluster model
      pred_dat = as.data.frame(mod_dat[,which(names(mod_dat)%in%clusters)])
      if (ncol(pred_dat)<2){
        colnames(pred_dat) = paste("`",clusters,"`",sep = "")
      }else{
        colnames(pred_dat) = paste("`",colnames(pred_dat),"`",sep = "")
      }
      #Predict the full dataset on the model
      plda <- predict(object = mod,
                      newdata = pred_dat)
      
      #Extract the position for the LDA for this group
      dataset = data.frame(Group = mod_dat[,"Cat"], lda = plda$x,Y = 0)
      
      #Get the group stats
      stats = dataset %>%
        group_by(Group) %>%
        summarise(Mean=mean(LD1), SD=sd(LD1))
      stats = as.data.frame(stats)
      stats$Model[3] = paste(mod[["xNames"]],collapse = "")
      stats$Naive_Healthy[3] = abs(stats[1,2] - stats[2,2])
      stats$Naive_Sick[3] = abs(stats[1,2] - stats[3,2])
      stats$Sick_Healthy[3] = abs(stats[3,2] - stats[2,2])
      
      #Add stats to the list
      table_stats[[index]] = stats
      
      #Draw density estimate
      ggplot(dataset, aes(x = LD1))+
        labs(y = paste("Denisity Estimate", sep=""))+
        xlim(-15,15)+
        geom_density(aes(fill = Group), alpha=0.75)+
        theme_bw()
      ggsave(paste(paste(clusters,collapse = " "),".jpeg",sep = ""))
      
      #Move the index up
      index = index +1
    }
  }
}

table_stats = do.call(rbind,table_stats)
write.xlsx(table_stats,paste("LDA_projection_stats.xlsx"))

