source("Spearman_Pearson_Comparison_function.R")

#Run the comparison
Run.Correlations(experiment1_file = "Grp 1 D35 Innate Cells 20181211 K=26.txt",
                 experiment2_file = "Grp 2 D35 Innate Cells 20181211 K=18.txt",
                 export_file = "DR3 D35 Innate cell Corr Comparison.csv")