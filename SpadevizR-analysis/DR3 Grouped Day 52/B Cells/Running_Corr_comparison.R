source("Spearman_Pearson_Comparison_function.R")

#Run the comparison
Run.Correlations(experiment1_file = "DR3 D52 Grp 1 B Cells 20190103 K=49.txt",
                 experiment2_file = "DR3 D52 Grp 2 B Cells 20190103 K=41.txt",
                 export_file = "DR3 D52 Bcell Corr Comparison.csv")