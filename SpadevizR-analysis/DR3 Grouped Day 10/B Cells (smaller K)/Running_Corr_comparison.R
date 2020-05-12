source("Spearman_Pearson_Comparison_function.R")

#Run the comparison
Run.Correlations(experiment1_file = "Grp 1 D10 B Cells 20181205 K=25.txt",
                 experiment2_file = "Grp 2 D11 B Cells 20181205 K=41.txt",
                 export_file = "DR3 D10 Bcell Corr Comparison.csv")