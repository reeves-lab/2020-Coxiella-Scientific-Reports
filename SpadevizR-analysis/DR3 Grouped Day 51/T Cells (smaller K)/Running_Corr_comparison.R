source("Spearman_Pearson_Comparison_function.R")

#Run the comparison
Run.Correlations(experiment1_file = "Grp 1 DR3 D51 T Cells 20181203 K = 29.txt",
                 experiment2_file = "Grp 2 DR3 D51 T Cells 20181203 K = 21.txt",
                 export_file = "DR3 D51 Tcell Corr Comparison.csv")