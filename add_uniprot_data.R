library(httr)
library(RCurl)
library(RJSONIO)
library(ggplot2)
library(Peptides)

rm(list=ls())

res_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/well_behaved/pumba_human_proteins_210525.txt")
res_table <- read.table(file=res_path, sep="\t", header = TRUE)

# internal library
source("./R/functions/all_functions.R")

nr_transmembrane_vec <- c()

# loop through table
for(i in 1:nrow(res_table)){
  print(paste0(i, " of ", nrow(res_table)))
  
  protein_ac <- res_table$protein.ac[i]
  uniprot_xml <- get_uniprot_xml(protein_ac)
  
  nr_transmembrane <- get_nr_transmembranes(uniprot_xml)
  nr_transmembrane_vec <- c(nr_transmembrane_vec, nr_transmembrane)
}

res_table <- cbind(res_table, nr_transmembrane_vec)
colnames(res_table) <- c(colnames(res_table), "transmembrane.nr")


  