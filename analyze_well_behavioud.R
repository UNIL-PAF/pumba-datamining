rm(list=ls())

# input
res_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/well_behaved/pumba_human_proteins_210518.txt")

# output
flt_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/well_behaved/pumba_human_proteins_firstAC_210518.txt")

res_table <- read.table(file=res_path, sep="\t", header = TRUE)
dim(res_table)


# keep only results where all samples got a result as the first ac
is_first_ac_cols <- grep("is.first.ac", colnames(res_table))
all_have_result <- apply(res_table[,is_first_ac_cols], 1, function(x) {if(any(x == FALSE) || any(is.na(x))) FALSE else TRUE })
res_table_2 <- res_table[all_have_result,]

write.table(res_table_2, file = flt_path, sep = "\t", row.names = FALSE)


# compare highest peak masses between replicates
mass_cols <- grep(".closest.peak.mass", colnames(res_table_2))

library("PerformanceAnalytics")
chart.Correlation(res_table_2[,mass_cols], histogram=TRUE, pch=19)
