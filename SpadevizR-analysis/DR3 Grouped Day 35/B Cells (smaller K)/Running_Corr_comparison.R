source("Spearman_Pearson_Comparison_function.R")

#Run the comparison
Run.Correlations(experiment1_file = "Grp 1 D35 B Cells 20181205 K=35.txt",
                 experiment2_file = "Grp 2 D35 B Cells 20181205 K=26.txt",
                 export_file = "DR3 D35 Bcell Corr Comparison.csv")