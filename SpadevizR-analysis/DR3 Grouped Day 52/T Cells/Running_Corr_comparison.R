source("Spearman_Pearson_Comparison_function.R")

#Run the comparison
Run.Correlations(experiment1_file = "DR3 D52 Grp 1 T Cells 20190103 K=33.txt",
                 experiment2_file = "DR3 D52 Grp 2 T Cells 20190103 K=25.txt",
                 export_file = "DR3 D52 Tcell Corr Comparison.csv")