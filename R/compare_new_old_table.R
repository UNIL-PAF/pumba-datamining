library(httr)
library(RCurl)
library(RJSONIO)
library(ggplot2)
library(Peptides)
library(diffdf)

rm(list=ls())

# new table
res_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/well_behaved/pumba_human_proteins_211118.txt")
res_table <- read.table(file=res_path, sep="\t", header = TRUE)


# load old table to compare
#old_res_path <- "/Users/rmylonas/tmp/datamining_pumba/results/well_behaved/pumba_human_proteins_210817.txt"
old_res_path <- "/Users/rmylonas/tmp/datamining_pumba/results/well_behaved/pumba_human_proteins_211005.txt"
old_res_table <- read.table(file=old_res_path, sep="\t", header = TRUE)

# transform all logical columns
old_colnames <- colnames(old_res_table)



for(k in 948:1000){
  
  print(k)
  
  # compare the results with the old table
  one_row <- res_table[k,]
  
  old_one_row <- old_res_table[k,]
  
  for(one_colname in old_colnames){
    if(is.na(one_row[one_colname]) && is.na(old_one_row[one_colname])) next
    diff <- diffdf( one_row[one_colname] , old_one_row[one_colname], suppress_warnings = TRUE, strict_numeric = FALSE)
    if(diffdf_has_issues(diff)){
      print("Difference in ")
      print(one_row$protein.ac)
      print(k)
      print(diff)
      #stop()
    }
  }
  
  
  
}

