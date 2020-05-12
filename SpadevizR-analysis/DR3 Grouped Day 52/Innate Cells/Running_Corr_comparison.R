source("Spearman_Pearson_Comparison_function.R")

#Run the comparison
Run.Correlations(experiment1_file = "DR3 D52 Grp 1 Innate Cells 20190103 K=41.txt",
                 experiment2_file = "DR3 D52 Grp 2 Innate Cells 20190103 K=33.txt",
                 export_file = "DR3 D52 Innate cell Corr Comparison.csv")