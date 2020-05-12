source("Spearman_Pearson_Comparison_function.R")

#Run the comparison
Run.Correlations(experiment1_file = "Grp 1 DR3 D51 Innate Cells 20180524 K = 35.txt",
                 experiment2_file = "Grp 2 DR3 D51 Innate Cells 20180524 K = 26.txt",
                 export_file = "DR3 D51 Innate cell Corr Comparison.csv")