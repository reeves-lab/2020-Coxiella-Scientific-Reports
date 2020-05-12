#Determine optimal module number for CNA - Qvax
#Joshua Hess
library(openxlsx)
library(factoextra)
library(NbClust)

#Read the correlation matrix
data_full = read.xlsx("D10-D51_Coefficients NCvVC Matrix.xlsx")

#Subset the data
data_clustering = data_full[1:228,230:238] #Indexes clusters as rows and clin measures as cols

#
results=NbClust(data_clustering,distance = "euclidean", min.nc = 5, max.nc = 20,
        method = "ward.D2")
results$Best.nc
results$Best.partition
print(results)
capture.output(results, file = "OptimalClusterSelection.txt")
