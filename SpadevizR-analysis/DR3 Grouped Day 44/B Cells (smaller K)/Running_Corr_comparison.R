source("Spearman_Pearson_Comparison_function.R")

#Run the comparison
Run.Correlations(experiment1_file = "20181205 Grp 1 DR3 D44 B Cells K = 43.txt",
                 experiment2_file = "20181205 Grp 2 DR3 D44 B Cells K = 26.txt",
                 export_file = "DR3 D44 Bcell Corr Comparison.csv")