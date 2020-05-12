# Classification-models
Scripts and functions for running multiple classification models (Feed forward neural network, linear discriminant analysis, support vector machine, naive bayes classifier) and exporting performance measures (classification accuracy with confidence intervals). Exhaustive search of parameters that maximize classification accuracy is performed and exported to excel file.

*Note: parameters (variables) are first filtered on the basis of absolute value of log odds coefficient from elastic net regression (see other folders)*

## Functions:
1) Read_data_ML
 - Function for reading and preparing data
2) ML_models_exhaustive
 - Function for running all models with exhaustive combinations

## Scripts:
1) Day10_ML_models
 - Script for Day 10 models
2) Day35_ML_models
 - Script for Day 35 models
3) Day44_ML_models
 - Script for Day 44 models
4) Day51_ML_models
 - Script for Day 51 models