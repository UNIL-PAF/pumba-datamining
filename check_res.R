library(httr)
library(RCurl)
library(RJSONIO)
library(ggplot2)
library(Peptides)


rm(list=ls())

res_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/well_behaved/pumba_human_proteins_210525.txt")
res_table <- read.table(file=res_path, sep="\t", header = TRUE)

# parameters
organism_param <- "human"
db_param <- "UP000005640_canonical_plus_splices_16092020"

# internal library
source("./R/functions/all_functions.R")
all_datasets <- get_all_datasets(organism=organism_param)

for(i in 1:length(all_datasets)){
 
  sample_id <- paste0(all_datasets[[i]]$sample, ".", all_datasets[[i]]$name)
  print(paste0("Missing ACs for ", sample_id, ": "))
  
  norm_path <- paste0("/Users/rmylonas/tmp/pumba/", all_datasets[[i]][["massFitResult"]][["proteinGroupsPath"]])  
  norm_table <- read.table(file=norm_path, sep="\t", header = TRUE) 
  norm_first_acs <- unlist(lapply(strsplit(as.character(norm_table$Majority.protein.IDs), ";"), function(x){x}))
  
  print(paste0("Nr rows: ", nrow(norm_table)))
  print(paste0("Nr acs: ", length(norm_first_acs)))
  
  first_ac <- which(res_table[,c(paste0(sample_id, ".nr.peaks"))] >= 1)
  res_ac <- res_table$protein.ac[first_ac]
  
  
  # what are the missing ones ?
  missing_acs <- setdiff(norm_first_acs, res_ac)
  
  # and in the other sense ?
  reverse_missing_acs <- setdiff(res_ac, norm_first_acs)
  
  print(paste0(length(missing_acs), " - ", length(reverse_missing_acs)))
}



# # HeLa 1160
# #dataset_id <- which(unlist(lapply(all_datasets, function(x){x$"name"})) == "11660")
# 
# 
# missing_row <- which(res_table$protein.ac == "A0A096LP01")
# res_table[missing_row,]
# 
# 
# # good, no ones are missing in this sense
# setdiff(res_ac_1160, norm_first_acs_1160)
# 
# 
# # look at the original data
# orig_path <- "/Users/rmylonas/Work/PAF/projects/pumba/data/maxquant/txt/proteinGroups.txt"
# orig_table <- read.table(file = orig_path, sep="\t", header = TRUE)
# orig_prot_row <- grep("A2VDF0", orig_table$Majority.protein.IDs)
# int_columns <- grep("Intensity.11660.", colnames(orig_table))
# orig_table[orig_prot_row,int_columns]
# 
# 
# 
# 
# colnames(res_table)
# 
# summary(res_table[,c("HeLa.11660.is.first.ac")])

# U2OS
u2os.path <- "/Users/rmylonas/Work/PAF/projects/pumba/data/maxquant/U2OS_9052_9053_9508_with_MBR.txt"
u2os <- read.table(file=u2os.path, sep="\t", header = TRUE)
nrow(u2os)
# 6897

# remove contaminants and rev
u2os.flt <- u2os[u2os$Potential.contaminant != "+", ]
# 6804
u2os.flt <- u2os[u2os$Reverse != "+", ]
# 6736

nr.peak.cols <- grep("U2OS.+.nr.peaks", colnames(res_table))
nr.peak.means <- apply(res_table[,nr.peak.cols],1,function(x){mean(x, na.rm=TRUE)})
sum(nr.peak.means > 0, na.rm = TRUE)
# 6414


sum(res_table$U2OS.9508.nr.peaks > 0, na.rm = TRUE)





