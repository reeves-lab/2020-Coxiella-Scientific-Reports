# Qvax-Mouse-Scientific-Reports
Analysis scripts for Qvax mouse manuscript submission to scientific reports. Link to paper:

## SpadevizR-analysis
Includes scripts for SpadevizR/cluster matching analysis of DR3 mice days 10-51. 

## Classification-models
Scripts and functions for running multiple classification models (Feed forward neural network, linear discriminant analysis, support vector machine, naive bayes classifier) and exporting performance measures (classification accuracy with confidence intervals). Exhaustive search of parameters that maximize classification accuracy is performed and exported to excel file.

*Note: parameters (variables) are first filtered on the basis of absolute value of log odds coefficient from elastic net regression (see other folders)*