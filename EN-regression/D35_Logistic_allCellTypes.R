

library(caret)
library(glmnet)
library(mlbench)
library(psych)
library(ggplot2)
library(openxlsx)
library(e1071) 
library(lattice)
library(fastDummies)
library(randomForest)
library(gridExtra)
library(dplyr)
library(stringr)
library(doParallel)

#----------------------Read Data - Change Sheet for different days--------------------------


dataD35T <- read.xlsx("HistoScore and Phenotype Data_T Cell Clusters(smaller K)_DR3.xlsx",
                      sheet = "Day 35")
dataD35T <- na.omit(dataD35T)
Category = dataD35T[c("Category")]
clusterD35T <- dplyr::select(dataD35T, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD35T)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD35T)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD35T)
prefix = "T D35 Cluster"
suffix = seq(1:n)
colnames(clusterD35T)[1:ncol(clusterD35T)] = paste(prefix, suffix, sep = " ")

dataD35B <- read.xlsx("HistoScore and Phenotype Data_B Cell Clusters(smallerK)_DR3.xlsx",
                      sheet = "Day 35")
dataD35B <- na.omit(dataD35B)
clusterD35B <- dplyr::select(dataD35B, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD35B)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD35B)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD35B)
prefix = "B D35 Cluster"
suffix = seq(1:n)
colnames(clusterD35B)[1:ncol(clusterD35B)] = paste(prefix, suffix, sep = " ")

dataD35I <- read.xlsx("HistoScore and Phenotype Data_Innate Cell Clusters(mixed K)_DR3.xlsx",
                      sheet = "Day 35")
dataD35I <- na.omit(dataD35I)
clusterD35I <- dplyr::select(dataD35I, starts_with("Cluster"))
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(clusterD35I)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(clusterD35I)[17:31] <- paste(prefix, suffix, sep = " ")

n = ncol(clusterD35I)
prefix = "I D35 Cluster"
suffix = seq(1:n)
colnames(clusterD35I)[1:ncol(clusterD35I)] = paste(prefix, suffix, sep = " ")


data = cbind(clusterD35T,
             clusterD35B,
             clusterD35I,
             Category)



#---------------------Data Rearrangement------------------------------------------------------

# data <- na.omit(data)
# # Category = data[c("Category")]
# data <- subset(data, select = -c(Total.Histoscore, Total.Antibody, Category))
# data
# rownames(data) = data$X
# data <- subset(data, select = -c(X))
# data = as.data.frame(lapply(data, as.numeric))        #Transform text to number
# data <- cbind(data, Category)
# str(data)
# #pairs.panels(data, lm = TRUE)


#Set up the dataframe for Binomial Lasso
#Predictors are scaled to mean = 0 and SD = 1
# cols <- c("HeartScore", "LiverScore", "SpleenScore", "LungScore")
# data[cols] <- lapply(data[cols], factor)  ## as.factor() could also be used

# dumsData <- dummy_cols(data, remove_first_dummy = FALSE)
# dumsData = dumsData[ ,-(ncol(dumsData))]
Clusters <- dplyr::select(data, contains("Cluster"))
# Continuous <- dplyr::select(dumsData, "Spleen.Burdon..copies.mg.", "Spleen.Body.wt...ratio.",
#                             "Day10.Antibody.OD", "Day24.Antibody.OD", "Day35.Antibody.OD")
# Ordinal <- dplyr::select(dumsData, "HeartScore_1", "LiverScore_2", "LiverScore_1", "SpleenScore_2",
#                          "SpleenScore_1", "LungScore_1")


#---------Both Continuous and Dummy Variables are scaled-----------------

# x <- data.matrix(data.frame(Clusters,Continuous, Ordinal), rownames.force = FALSE)
x <- data.matrix(data.frame(Clusters), rownames.force = FALSE)
x <- scale(x, center = TRUE, scale = TRUE)
x = as.data.frame(x)
prefix = "Grp1 Mouse"
suffix = seq(1,16)
rownames(x)[1:16] <- paste(prefix, suffix, sep = " ")
prefix = "Grp2 Mouse"
suffix = seq(2,16)
rownames(x)[17:31] <- paste(prefix, suffix, sep = " ")
# 
# n = ncol(Clusters)
# prefix = "Cluster"
# suffix = seq(1:n)
# colnames(x)[1:ncol(Clusters)] = paste(prefix, suffix, sep = " ")
# colnames(x)[(ncol(Clusters)+1):(ncol(Clusters)+5)] = c("Spleen Burden (copies/mg)", 
#                                                        "Spleen/Body wt. (ratio)",
#                                                        "Day10 Antibody OD",
#                                                        "Day24 Antibody OD",
#                                                        "Day35 Antibody OD")
# colnames(x)[(ncol(Clusters)+6):(ncol(Clusters)+11)] = c("Heart score (high)",
#                                                         "Liver score (high)",
#                                                         "Liver score (medium)",
#                                                         "Spleen score (high)",
#                                                         "Spleen score (medium)",
#                                                         "Lung score (high)")
x = cbind(x, Category)

# NU = x[x$Category == "Naive Unchallenge", ]
# NC = x[x$Category == "Naive Challenge", ]
# VU = x[x$Category == "Vaccinated Unchallenge", ]
# VC = x[x$Category == "Vaccinated Challenge", ]


registerDoParallel(cores = 4)

##############################
# Naive Unchal vs Naive Chal
##############################


data = x
data$Category = as.factor(data$Category)
data$Category
levels(data$Category) <- make.names(levels(factor(data$Category)))



#--------------------Custom Control Parameters---------------------------

custom <- caret::trainControl(method = "repeatedcv",
                              number = 4,
                              repeats = 100,
                              verboseIter = TRUE, classProbs = TRUE)



#-------------------------LASSO-----------------------------------------

set.seed(433)
lasso <- caret::train(Category ~ . ,
                      data,
                      method = 'glmnet', family = "binomial", standardize = FALSE,
                      tuneGrid = expand.grid(
                        lambda = seq(0.00000001, 1,length = 100),
                        alpha = 1),
                      trControl = custom)



#-----------Plot Cross-Validation ------------------------------

plot1 = ggplot(data = lasso) +theme(text = element_text(size=17),              #X and Y axis title sizes 
                                    axis.text.x = element_text(face="bold", color="black", size=16, angle=0),    #X tick Size
                                    axis.text.y = element_text(face="bold", color="black", size=16, angle=0),    #Y tick Size
                                    panel.background = element_blank(),                                          #Remove grey background
                                    axis.line = element_line(colour = "black"),                                  #Black line at start of axis
                                    panel.grid = element_line(colour = "grey"),                                  #Grid lines
                                    plot.margin = unit(c(2,2,2,2), "cm"),                                        #Plot margins
                                    axis.title.y = element_text(angle = 90, vjust = 10, hjust = 0.5),            #Y axis title position
                                    axis.title.x = element_text(angle = 0, vjust = -10, hjust = 0.5))            #X axis title position
ggsave("LogLasso_CV_NvV.tiff")



#-----------Plot Importance ------------------------------

w = varImp(lasso, scale = FALSE)
w = w$importance
a = rownames(w)
a = str_replace_all(a, "[.]", " ")
rownames(w) = a

w = w[order(-w$Overall), , drop = FALSE]
values = w$Overall
LABEL = rownames(w)
LABEL = factor(LABEL, levels = LABEL[order(values)])
plot2 = ggplot(w) + 
  geom_bar(stat='identity', aes(fill=values, x = LABEL, y = values), fill = "black", width = 0.5) + theme(text = element_text(size=17),                         
                                                                                                          axis.text.x = element_text(face="bold", color="black", size=1.8, angle=0),
                                                                                                          axis.text.y = element_text(face="bold", color="black", size=1.8, angle=0),
                                                                                                          panel.background = element_blank(),
                                                                                                          axis.line = element_line(colour = "black", size = 0.15),
                                                                                                          panel.grid = element_line(colour = "grey", size = 0.1),
                                                                                                          plot.margin = unit(c(2.5,2.5,2.5,2.5), "mm"),
                                                                                                          axis.title.y = element_text(angle = 90, vjust = 0.4, hjust = 0.5, size = 4),
                                                                                                          axis.title.x = element_text(angle = 0, vjust = -0.25, hjust = 0.5, size = 4)) +
  labs(x = "Feature", y = "Importance (|logit|)")+
  coord_flip()
ggsave("LogLasso_Importance_NvV.tiff", plot = plot2, units = "mm", width = 40, height = 55)


#-----------Plot Coefficient Plots ------------------------------
tiff("LogLasso_Coefficients_Lambda_NvV.tiff", width = 450, height = 350)
plot(lasso$finalModel, xvar = 'lambda', label = TRUE)
dev.off()

tiff("LogLasso_Coefficients_Deviance_NvV.tiff", width = 450, height = 350)
plot(lasso$finalModel, xvar = "dev", label = TRUE )
dev.off()




#-----------Export Data Tables ------------------------------

cvSummary = as.data.frame(lasso$results)
write.xlsx(cvSummary, "cvSummary_Log Lasso_NvV.xlsx")



best <- lasso$finalModel
Coefficients = as.matrix(coef(best, s = lasso$bestTune$lambda))
b = rownames(Coefficients)
b = str_replace_all(b, "[.]", " ")
rownames(Coefficients) = b
Coefficients = as.data.frame(Coefficients)
colnames(Coefficients)[1] = "Change in Log Odds per Unit Increase"

# n = ncol(Clusters)
# prefix = "Cluster"
# suffix = seq(1:n)
# rownames(Coefficients)[2:(ncol(Clusters)+1)] = paste(prefix, suffix, sep = " ")
# rownames(Coefficients)[(ncol(Clusters)+2):(ncol(Clusters)+6)] = c("Spleen Burden (copies/mg)", 
#                                                                   "Spleen/Body wt. (ratio)",
#                                                                   "Day10 Antibody OD",
#                                                                   "Day24 Antibody OD",
#                                                                   "Day35 Antibody OD")
# rownames(Coefficients)[(ncol(Clusters)+7):(ncol(Clusters)+12)] = c("Heart score (high)",
#                                                                    "Liver score (high)",
#                                                                    "Liver score (medium)",
#                                                                    "Spleen score (high)",
#                                                                    "Spleen score (medium)",
#                                                                    "Lung score (high)")
Coefficients
write.xlsx(Coefficients, "Coefficients_LogLasso_LogOdds_NvV.xlsx", row.names = TRUE)



#calculating odds ratios
OddsRatio = as.matrix(exp(Coefficients))
OddsRatio = as.data.frame(OddsRatio)
colnames(OddsRatio)[1] = "Odds ratio between Unit 0 and Unit 1 (Unit Increase)"
OddsRatio
write.xlsx(OddsRatio, "Coefficients_LogLasso_OddsRatios_NvV.xlsx", row.names = TRUE)

#calculating probabilites
prob <- exp(Coefficients)/(1+exp(Coefficients))
Probabilities = as.data.frame(prob)
colnames(Probabilities)[1] <- "Probability associated with a one Unit Increase (To get Change, use 0.5 as reference)"
Probabilities
write.xlsx(Probabilities, "Coefficients_LogLasso_Probabilities_NvV.xlsx", row.names = TRUE)


#calculating odds ratios
Coefficients
exp(Coefficients)

#calculating probabilites
prob <- exp(Coefficients)/(1+exp(Coefficients))
prob
a = as.data.frame(prob)
colnames(a)[1] <- "Probability"



#----------------Plotting Log Odds Plot-------------------------------------------------

Coefficients = subset(Coefficients, Coefficients$`Change in Log Odds per Unit Increase`!=0)
Coefficients = Coefficients[2:nrow(Coefficients),,drop = FALSE]
Coefficients = Coefficients[order(-Coefficients$`Change in Log Odds per Unit Increase`), , drop = FALSE]

Coefficients$`Group` <- rownames(Coefficients)  # create new column for car names
Coefficients$Group <- ifelse(Coefficients$`Change in Log Odds per Unit Increase` < 0, "Naive", "Vaccinated")  # above / below avg flag
Coefficients

# Diverging Barchart
logoddsvalue = Coefficients$`Change in Log Odds per Unit Increase`
features = rownames(Coefficients)
features = factor(features, levels = features[order(logoddsvalue)])
axval = (max(abs(logoddsvalue)) + 0.1)

plot5 = ggplot(Coefficients, aes(x=features, y=logoddsvalue)) + 
  geom_bar(stat='identity', aes(fill=Group), width=.4)  +
  scale_fill_manual(name="Group", 
                    labels = c("Naive", "Vaccinated"), 
                    values = c("Naive"="black", "Vaccinated"="grey60"), drop = FALSE) + 
  labs(subtitle="", 
       title= "") + 
  theme(text = element_text(size=0.5),                         
        axis.text.x = element_text(face="bold", color="black", size=2.5, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=2.5, angle=0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.15),
        panel.grid = element_line(colour = "grey", size = 0.1),
        plot.margin = unit(c(2,1,1,2), "mm"),
        axis.title.y = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 4),
        axis.title.x = element_text(angle = 0, vjust = -0.3, hjust = 0.5, size = 4),
        plot.title = element_text(vjust=10, hjust = 0),
        plot.subtitle = element_text(vjust = 10),
        legend.title = element_text(size=2), 
        legend.text = element_text(size=2),
        legend.key.size = unit(0.3,"line"),
        legend.margin = ggplot2::margin(t=0.5, unit = "mm")) +
  labs(x = "Feature", y = expression(paste(Delta, "Log Odds/","Unit Increase"))) + 
  scale_y_continuous(limits = c(-axval, axval)) +
  coord_flip()
ggsave("LogLasso_LogOddsPlot_NvV.tiff", plot = plot5, units = "mm", width = 40, height = 35)




#----------------Plotting Probability Plot-------------------------------------------------

for (i in 1:nrow(a)) {
  if(a$Probability[i] < 0.5){
    a$Probability1[i] <- (-(1 - a$Probability[i]))
  }else{
    a$Probability1[i] <- a$Probability[i]
  }
}

a = as.data.frame(a[,-1])
rownames_data_a <- rownames(prob)
rownames(a) <- rownames_data_a
colnames(a)[1] <- "Probability"

for (i in 1:nrow(a)){
  if(a$Probability[i] > 0){
    a$Probability[i] = (a$Probability[i] - 0.5)
  }else{
    a$Probability[i] = (a$Probability[i] + 0.5)
  }
}


a = subset(a, Probability!=0)
a = a[2:nrow(a),,drop = FALSE]
a = a[order(-a$Probability), , drop = FALSE]

a$`Group` <- rownames(a)  # create new column for car names
a$Group <- ifelse(a$Probability < 0, "Naive", "Vaccinated")  # above / below avg flag
a

# Diverging Barchart
probvalues = a$Probability
feat = rownames(a)
feat = factor(feat, levels = feat[order(probvalues)])
axval = (0.5)


plot6 = ggplot(a, aes(x=feat, y=probvalues)) + 
  geom_bar(stat='identity', aes(fill=Group), width=.4)  +
  scale_fill_manual(name="Group", 
                    labels = c("Naive", "Vaccinated"), 
                    values = c("Vaccinated"="grey60", "Naive"="black")) + 
  labs(subtitle="", 
       title= "") + 
  theme(text = element_text(size=0.5),                         
        axis.text.x = element_text(face="bold", color="black", size=2.5, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=2.5, angle=0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid = element_line(colour = "grey", size = 0.1),
        plot.margin = unit(c(2,1,1,2), "mm"),
        axis.title.y = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 4),
        axis.title.x = element_text(angle = 0, vjust = -0.3, hjust = 0.5, size = 4),
        plot.title = element_text(vjust=10, hjust = 0),
        plot.subtitle = element_text(vjust = 10),
        legend.title = element_text(size=2), 
        legend.text = element_text(size=2),
        legend.key.size = unit(0.3,"line"),
        legend.margin = ggplot2::margin(t=0.5, unit = "mm")) +
  labs(x = "Feature", y = expression(paste(Delta, "Probability/Unit Increase (0 to 1)"))) +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  coord_flip()
ggsave("LogLasso_ProbabilityPlot_NvV.tiff", plot = plot6, units = "mm", width = 40, height = 35)





# #----------------Test model on unseen data (Test data) ------------------
# 
# pred = predict(lasso, newdata=test)
# accuracy <- table(pred, test[,"Category"])
# sum(diag(accuracy))/sum(accuracy)
# AccuracyTest = sum(diag(accuracy))/sum(accuracy)
# write.xlsx(AccuracyTest, "Log Lasso_AccuracyTest_NvV_ClustersHisto.xlsx")

#-------------------------------------------------------------------------






#--------------------Custom Control Parameters---------------------------

custom <- caret::trainControl(method = "repeatedcv",
                              number = 4,
                              repeats = 100,
                              verboseIter = TRUE, classProbs = TRUE)




#----------------Elastic Net----------------------------------------------

set.seed(433)
ElasticNet <- caret::train(Category ~ . ,
                           data,
                           method = 'glmnet', family = "binomial", standardize = FALSE,
                           tuneGrid = expand.grid(
                             lambda = seq(0.00000001, 1,length = 100),
                             alpha = seq(0,1,length = 50)),
                           trControl = custom)



#-----------Plot Cross-Validation ------------------------------

plot1 = ggplot(data = ElasticNet) +theme(text = element_text(size=17),              #X and Y axis title sizes 
                                         axis.text.x = element_text(face="bold", color="black", size=23, angle=0),    #X tick Size
                                         axis.text.y = element_text(face="bold", color="black", size=23, angle=0),    #Y tick Size
                                         panel.background = element_blank(),                                          #Remove grey background
                                         axis.line = element_line(colour = "black"),                                  #Black line at start of axis
                                         panel.grid = element_line(colour = "grey"),                                  #Grid lines
                                         plot.margin = unit(c(2,2,2,2), "cm"),                                        #Plot margins
                                         axis.title.y = element_text(angle = 90, vjust = 10, hjust = 0.5, size = 25),            #Y axis title position
                                         axis.title.x = element_text(angle = 0, vjust = -10, hjust = 0.5, size = 25),       #X axis title position
                                         legend.title = element_text(size=19), 
                                         legend.text = element_text(size=19))
ggsave("LambdaCV_EN_NvV.tiff", plot = plot1, units = "mm", width = 500, height = 300, limitsize = FALSE)




#-----------Plot Importance ------------------------------

w = varImp(ElasticNet, scale = FALSE)
w = w$importance
a = rownames(w)
a = str_replace_all(a, "[.]", " ")
rownames(w) = a

w = w[order(-w$Overall), , drop = FALSE]
values = w$Overall
LABEL = rownames(w)
LABEL = factor(LABEL, levels = LABEL[order(values)])
plot2 = ggplot(w) + 
  geom_bar(stat='identity', aes(fill=values, x = LABEL, y = values), fill = "black", width = 0.5) + theme(text = element_text(size=17),                         
                                                                                                          axis.text.x = element_text(face="bold", color="black", size=1.8, angle=0),
                                                                                                          axis.text.y = element_text(face="bold", color="black", size=1.8, angle=0),
                                                                                                          panel.background = element_blank(),
                                                                                                          axis.line = element_line(colour = "black", size = 0.15),
                                                                                                          panel.grid = element_line(colour = "grey", size = 0.1),
                                                                                                          plot.margin = unit(c(2.5,2.5,2.5,2.5), "mm"),
                                                                                                          axis.title.y = element_text(angle = 90, vjust = 0.4, hjust = 0.5, size = 4),
                                                                                                          axis.title.x = element_text(angle = 0, vjust = -0.25, hjust = 0.5, size = 4)) +
  labs(x = "Feature", y = "Importance (|logit|)")+
  coord_flip()
ggsave("LogEN_Importance_NvV.tiff", plot = plot2, units = "mm", width = 40, height = 55)


#-----------Plot Coefficient Plots ------------------------------
tiff("LogEN_Coefficients_Lambda_NvV.tiff", width = 450, height = 350)
plot(ElasticNet$finalModel, xvar = 'lambda', label = TRUE)
dev.off()

tiff("LogEN_Coefficients_Deviance_NvV.tiff", width = 450, height = 350)
plot(ElasticNet$finalModel, xvar = "dev", label = TRUE )
dev.off()



#-----------Export Data Tables ------------------------------

cvSummaryEN = as.data.frame(ElasticNet$results)
write.xlsx(cvSummaryEN, "cvSummary_LogElasticNet_NvV.xlsx")



bestEN <- ElasticNet$finalModel
CoefficientsEN = as.matrix(coef(bestEN, s = ElasticNet$bestTune$lambda))
c = rownames(CoefficientsEN)
c = str_replace_all(c, "[.]", " ")
rownames(CoefficientsEN) = c
CoefficientsEN = as.data.frame(CoefficientsEN)
colnames(CoefficientsEN)[1] = "Change in Log Odds per Unit Increase"
# n = ncol(Clusters)
# prefix = "Cluster"
# suffix = seq(1:n)
# rownames(CoefficientsEN)[2:(ncol(Clusters)+1)] = paste(prefix, suffix, sep = " ")
# rownames(CoefficientsEN)[(ncol(Clusters)+2):(ncol(Clusters)+6)] = c("Spleen Burden (copies/mg)", 
#                                                                     "Spleen/Body wt. (ratio)",
#                                                                     "Day10 Antibody OD",
#                                                                     "Day24 Antibody OD",
#                                                                     "Day35 Antibody OD")
# rownames(CoefficientsEN)[(ncol(Clusters)+7):(ncol(Clusters)+12)] = c("Heart score (high)",
#                                                                      "Liver score (high)",
#                                                                      "Liver score (medium)",
#                                                                      "Spleen score (high)",
#                                                                      "Spleen score (medium)",
#                                                                      "Lung score (high)")
CoefficientsEN
write.xlsx(CoefficientsEN, "Coefficients_LogEN_LogOdds_NvV.xlsx", row.names = TRUE)


#calculating odds ratios
OddsRatioEN = as.matrix(exp(CoefficientsEN))
OddsRatioEN = as.data.frame(OddsRatioEN)
colnames(OddsRatioEN)[1] = "Odds ratio between Unit 0 and Unit 1 (Unit Increase)"
OddsRatioEN
write.xlsx(OddsRatioEN, "Coefficients_LogEN_OddsRatios_NvV.xlsx", row.names = TRUE)


#calculating probabilites
probEN <- exp(CoefficientsEN)/(1+exp(CoefficientsEN))
ProbabilitiesEN = as.data.frame(probEN)
colnames(ProbabilitiesEN)[1] <- "Probability associated with a one Unit Increase (To get Change, use 0.5 as reference)"
ProbabilitiesEN
write.xlsx(ProbabilitiesEN, "Coefficients_LogEN_Probabilities_NvV.xlsx", row.names = TRUE)


#calculating odds ratios
CoefficientsEN
exp(CoefficientsEN)

#calculating probabilites
probEN <- exp(CoefficientsEN)/(1+exp(CoefficientsEN))
probEN
a = as.data.frame(probEN)
colnames(a)[1] <- "Probability"



#----------------Plotting Log Odds Plot-------------------------------------------------

CoefficientsEN = subset(CoefficientsEN, CoefficientsEN$`Change in Log Odds per Unit Increase`!=0)
CoefficientsEN = CoefficientsEN[2:nrow(CoefficientsEN),,drop = FALSE]
CoefficientsEN = CoefficientsEN[order(-CoefficientsEN$`Change in Log Odds per Unit Increase`), , drop = FALSE]

CoefficientsEN$`Group` <- rownames(CoefficientsEN)  # create new column for car names
CoefficientsEN$Group <- ifelse(CoefficientsEN$`Change in Log Odds per Unit Increase` < 0, "Naive", "Vaccinated")  # above / below avg flag
CoefficientsEN


# Diverging Barchart
logoddsvalue = CoefficientsEN$`Change in Log Odds per Unit Increase`
features = rownames(CoefficientsEN)
features = factor(features, levels = features[order(logoddsvalue)])
axval = (max(abs(logoddsvalue)) + 0.1)

plot5 = ggplot(CoefficientsEN, aes(x=features, y=logoddsvalue)) + 
  geom_bar(stat='identity', aes(fill=Group), width=.6)  +
  scale_fill_manual(name="Group", 
                    labels = c("Naive", "Vaccinated"), 
                    values = c("Vaccinated"="grey60", "Naive"="black")) + 
  labs(subtitle="", 
       title= "") + 
  theme(text = element_text(size=0.5),                         
        axis.text.x = element_text(face="bold", color="black", size=2, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=2, angle=0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.15),
        panel.grid = element_line(colour = "grey", size = 0.1),
        plot.margin = unit(c(2,1,1,2), "mm"),
        axis.title.y = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 4),
        axis.title.x = element_text(angle = 0, vjust = -0.3, hjust = 0.5, size = 4),
        plot.title = element_text(vjust=10, hjust = 0),
        plot.subtitle = element_text(vjust = 10),
        legend.title = element_text(size=2), 
        legend.text = element_text(size=2),
        legend.key.size = unit(0.3,"line"),
        legend.margin = ggplot2::margin(t=-0.5, unit = "mm")) +
  labs(x = "Feature", y = expression(paste(Delta, "Log Odds/","Unit Increase"))) + 
  scale_y_continuous(limits = c(-axval, axval)) +
  coord_flip()
ggsave("LogEN_LogOddsPlot_NvV.tiff", plot = plot5, units = "mm", width = 55, height = 50)




#-------------------------------------------------------------------------------------

#Plotting the results for Probabilities Elastic Net
for (i in 1:nrow(a)) {
  if(a$Probability[i] < 0.5){
    a$Probability1[i] <- (-(1 - a$Probability[i]))
  }else{
    a$Probability1[i] <- a$Probability[i]
  }
}

a = as.data.frame(a[,-1])
rownames_data_a <- rownames(prob)
rownames(a) <- rownames_data_a
colnames(a)[1] <- "Probability"

for (i in 1:nrow(a)){
  if(a$Probability[i] > 0){
    a$Probability[i] = (a$Probability[i] - 0.5)
  }else{
    a$Probability[i] = (a$Probability[i] + 0.5)
  }
}



a = subset(a, Probability!=0)
a = a[2:nrow(a),,drop = FALSE]
a = a[order(-a$Probability), , drop = FALSE]

a$`Group` <- rownames(a)  # create new column for car names
a$Group <- ifelse(a$Probability < 0, "Naive", "Vaccinated")  # above / below avg flag
a

# Diverging Barchart
probvalues = a$Probability
feat = rownames(a)
feat = factor(feat, levels = feat[order(probvalues)])
axval = (0.5)

plot6 = ggplot(a, aes(x=feat, y=probvalues)) + 
  geom_bar(stat='identity', aes(fill=Group), width=.4)  +
  scale_fill_manual(name="Group", 
                    labels = c("Naive", "Vaccinated"), 
                    values = c("Vaccinated"="grey60", "Naive"="black")) + 
  labs(subtitle="", 
       title= "") + 
  theme(text = element_text(size=0.5),                         
        axis.text.x = element_text(face="bold", color="black", size=2, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=2, angle=0),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.15),
        panel.grid = element_line(colour = "grey", size = 0.1),
        plot.margin = unit(c(2,1,1,2), "mm"),
        axis.title.y = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 4),
        axis.title.x = element_text(angle = 0, vjust = -0.3, hjust = 0.5, size = 4),
        plot.title = element_text(vjust=10, hjust = 0),
        plot.subtitle = element_text(vjust = 10),
        legend.title = element_text(size=2), 
        legend.text = element_text(size=2),
        legend.key.size = unit(0.3,"line"),
        legend.margin = ggplot2::margin(t=0.5, unit = "mm")) +
  labs(x = "Feature", y = expression(paste(Delta, "Probability/Unit Increase (0 to 1)"))) +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  coord_flip()
ggsave("LogEN_ProbabilityPlot_NvV.tiff", plot = plot6, units = "mm", width = 55, height = 50)



# #----------------Test Elastic Net model on unseen data (Test data) ------------------
# 
# predEN = predict(ElasticNet, newdata=test)
# accuracyEN <- table(predEN, test[,"Category"])
# sum(diag(accuracyEN))/sum(accuracyEN)
# AccuracyTestEN = sum(diag(accuracyEN))/sum(accuracyEN)
# write.xlsx(AccuracyTestEN, "LogElasticNet_AccuracyTest_NvV_ClustersHisto.xlsx")



#-------------Compare Lasso and Elastic Net--------------------------------------------

model_list <- list(EN = ElasticNet, LASSO = lasso)
res <- resamples(model_list)
summary(res)
SUM = as.data.frame(summary(res)$statistics)
write.xlsx(res$values, "LogENvLasso_ResamplesFull_NvV.xlsx", row.names=TRUE)
write.xlsx(SUM, "LogENvLasso_ResampleSummary_NvV.xlsx", row.names=TRUE)

SUM = as.matrix(SUM)
tiff("LogENvLasso_ResampleSummary_NvV.tiff", width = 2000, height = 800)
barplot(SUM, beside = TRUE, legend.text = TRUE, cex.axis = 1.5, cex.names = 1.25)
dev.off()













