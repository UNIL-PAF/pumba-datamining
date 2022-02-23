library(httr)
library(RCurl)
library(RJSONIO)
library(ggplot2)
library(Peptides)


rm(list=ls())

# results
res_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/slice_mass_limits/human_slice_mass_limits_210930.txt")

# parameters
organism_param <- "human"
db_param <- "UP000005640_canonical_plus_splices_16092020"

source("./R/functions/all_functions.R")

all_datasets <- get_all_datasets(organism=organism_param)

my_colnames <- c("cell_line", "replicate", "slice", "start_mass", "end_mass", "slice_mass")
res_table <- data.frame(matrix(ncol=length(my_colnames),nrow=0, dimnames=list(NULL, my_colnames)))

# loop through datasets
for(k in 1:length(all_datasets)){
  one_dataset <- all_datasets[[k]]
  
  cell_line <- one_dataset$sample
  replicate <- one_dataset$name

  mass_fit_params <- one_dataset$massFitResult$massFitCoeffs
  fit_func <- function(x) { mass_fit_params[1] + mass_fit_params[2] * x + mass_fit_params[3] * x^2 + mass_fit_params[4] * x^3 }
  
  # loop through slides
  for(i in 1:length(one_dataset$massFitResult$massFits)){
     start_mass <- fit_func(i - 0.5)
     end_mass <- fit_func(i + 0.5)
     slice_mass <- fit_func(i)
     slice <- i
     
     one_row <- c(cell_line, replicate, slice, 10^end_mass, 10^start_mass, 10^slice_mass)
     res_table <- rbind(res_table, one_row)
  }
}

colnames(res_table) <- my_colnames

write.table(res_table, file = res_path, sep = "\t", row.names = FALSE)
