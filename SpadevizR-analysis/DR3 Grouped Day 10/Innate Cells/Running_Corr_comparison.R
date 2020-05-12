source("Spearman_Pearson_Comparison_function.R")

#Run the comparison
Run.Correlations(experiment1_file = "Grp 1 D10 Innate Cells 20180706 K=25.txt",
                 experiment2_file = "Grp 2 D11 Innate Cells 20180706 K=33.txt",
                 export_file = "DR3 D10 Innate cell Corr Comparison.csv")