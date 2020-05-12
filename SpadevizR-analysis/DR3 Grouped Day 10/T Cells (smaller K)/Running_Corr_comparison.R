source("Spearman_Pearson_Comparison_function.R")

#Run the comparison
Run.Correlations(experiment1_file = "Grp 1 D10 T Cells 20181204 K=17.txt",
                 experiment2_file = "Grp 2 D10 T Cells 20181204 K=25.txt",
                 export_file = "DR3 D10 Tcell Corr Comparison.csv")