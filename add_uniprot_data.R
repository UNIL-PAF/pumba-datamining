library(httr)
library(RCurl)
library(RJSONIO)
library(ggplot2)
library(Peptides)

rm(list=ls())

output_path <- "/Users/rmylonas/tmp/datamining_pumba/results/well_behaved/pumba_human_proteins_211005_PTMs.txt"

res_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/well_behaved/pumba_human_proteins_211005.txt")
res_table <- read.table(file=res_path, sep="\t", header = TRUE)

# internal library
source("./R/functions/all_functions.R")

ptms <- c()
transmembranes <- c()
glycosylations <- c()
compositionally_biased_regions <- c()
ubiquitins <- c()

# loop through table
for(i in 1:nrow(res_table)){
  print(paste0(i, " of ", nrow(res_table)))
  
  protein_ac <- res_table$protein.ac[i]
  uniprot_xml <- get_uniprot_xml(protein_ac)
  
  transmembrane <- get_nr_transmembranes(uniprot_xml)
  transmembranes <- c(transmembranes, transmembrane)
  
  glycosylation <- get_glycosylations(uniprot_xml)
  glycosylations <- c(glycosylations, glycosylation)
  
  compositionally_biased_region <- get_compositionally_biased_region(uniprot_xml)
  compositionally_biased_regions <- c(compositionally_biased_regions, compositionally_biased_region)
  
  ptm <- get_ptms(uniprot_xml)
  ptms <- c(ptms, ptm)  
  
  ubiquitin <- get_ubiquitin(uniprot_xml, i)
  ubiquitins <- c(ubiquitins, ubiquitin)
  
  #get_is_cleaved(uniprot_xml)
}

res_table <- cbind(res_table, ptms, transmembranes, glycosylations, compositionally_biased_regions, ubiquitins)
colnames(res_table) <- c(head(colnames(res_table), -5), "ptms", "transmembrane", "glycosylations", "compositionally.biased.region", "ubiquitin")

write.table(res_table, file = output_path, sep = "\t", row.names = FALSE)
  