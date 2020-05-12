# Qvax-Mouse-Scientific-Reports
Analysis scripts for Qvax mouse manuscript submission to scientific reports. Link to paper:

## SpadevizR-analysis
Includes scripts for SpadevizR/cluster matching analysis of DR3 mice days 10-52. 

## Classification-models
Scripts and functions for running multiple classification models (Feed forward neural network, linear discriminant analysis, support vector machine, naive bayes classifier) and exporting performance measures (classification accuracy with confidence intervals). Exhaustive search of parameters that maximize classification accuracy is performed and exported to excel file.

*Note: parameters (variables) are first filtered on the basis of absolute value of log odds coefficient from elastic net regression (see other folders)*

## Phenotypic-classes
Analysis scripts for phenotypic classes derived from hierarchical clustering and heatmap representation.

## LDA
Linear discriminant analysis per day (folder) of DR3 mice. Analysis considered combinations of 3 clusters for predicting class assignment and subsequent projection of results as a density distribution.

*Note: parameters (variables) are first filtered on the basis of absolute value of log odds coefficient from elastic net regression (see other folders)*

## EN-regression
Elastic net (and lasso comparison) per day (10-52) of all cell types for predicting class assignment of mice. Waterfall plots included, along with cross validation summary plotting.

## Correlation-analysis
Cluster frequency to clinical measure correlation analysis with module detection, heatmap generation, and t-SNE creation. T-sne plots included that map elastic net clusters to t-SNE plot.
