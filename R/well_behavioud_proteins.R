library(httr)
library(ggplot2)
library(Peptides)

rm(list=ls())

# results
res_path <- ("/Users/admin/tmp/datamining_pumba/results/well_behaved/well_behaved_prots.txt")

# threshold pour l'intensité totale pour cette protein
good_peak_threshold <- 0.8

# internal library
source("./R/functions/all_functions.R")

all_datasets <- get_all_datasets()

#nr_proteins <- 100
#plot_proteins <- c() #"H7C417", "O60307", "E9PLN8", "E7EVH7", "A2RTX5", "A2RRP1", "A0A0A6YYH1")

all_protein_groups <- get_all_protein_groups(all_datasets)

# just take the protein_acs from the first dataset
protein_acs <- get_protein_acs(all_protein_groups[[1]])

# loop over all proteins
nr_prot_loop <- if(! exists("nr_proteins")) length(protein_acs) else nr_proteins

dataset_ids <- get_dataset_ids(all_datasets)

all_protein_groups <- get_all_protein_groups(all_datasets)

res_table <- data.frame(stringsAsFactors=FALSE)
valid_slices_array <- array()
valid_entries <- 0
invalid_entries <- 0

for(k in 1:nr_prot_loop){
  
  print(paste0(k, " of ", nr_prot_loop))
  protein_ac <- protein_acs[k]
  
  # theo mol weight
  protein <- all_protein_groups[[1]][k,]
  theo_mol_weight <- as.numeric(protein$Mol..weight..kDa.)
  
  #k <- which(protein_acs == "P02786")
  row <- data.frame(protein_ac, theo_mol_weight)
  
  valid_entry <- TRUE
  
  for(i in 1:length(dataset_ids)){
  
    dataset_id <- dataset_ids[i]
    
    # get merged data from backend or cache
    protein_merges <- get_single_protein_merge(protein_ac, dataset_id)
    
    # for some datasets the protein is not detected
    if(length(protein_merges) == 0){
      valid_entry <- FALSE
      break
    }
    
    # since we load only one dataset at the time there is only 1 list
    protein_merge <- protein_merges[[1]]
    
    mass <- get_masses(protein_merge)
    ints <- get_ints(protein_merge)
    
    # get the indexes of intensities corresponding to a peak (see "config.R" for thresholds)
    peak_idxs <- get_peak_indexes(ints)
    
    # just skip with a warning if no peak was detected
    if(length(peak_idxs) < 1){
      valid_entry <- FALSE
      #warning(paste0("Could not find any peak in [", protein_ac, "]."))
      break
    }
    
    peaks_masses <- mass[peak_idxs]
    peaks_ints <- ints[peak_idxs]
    
    highest_peak <- which(peaks_ints == max(peaks_ints, na.rm=TRUE))
    highest_peak_idx <- peak_idxs[highest_peak]
    
    peak_limits <- get_peak_limits(ints, highest_peak_idx)
    peak_limits_masses <- mass[peak_limits]
  
    # find the correct limits in the slices
    log_mass_fits <- all_datasets[[i]]$massFitResult$massFits  
    mass_fits <- 10^(unlist(log_mass_fits))
    
    # find the slices which lie within the limits
    valid_slices <- which(mass_fits >= peak_limits_masses[1] & mass_fits <= peak_limits_masses[2])
    
    # check if this makes 80% of the total intensity
    slice_ints <- unlist(protein_merge$proteins[[1]]$intensities)
    slice_ints_sum <- sum(slice_ints[valid_slices])
    good_peak <- (slice_ints_sum / sum(slice_ints)) > good_peak_threshold
    if(! good_peak){
      valid_entry <- FALSE
      #warning(paste0("Bad peak for protein [", protein_ac, "] in dataset [", paste(all_datasets[[i]]$sample, all_datasets[[i]]$name, sep=".")  ,"]."))
      #warning(paste0("percentage was [", (slice_ints_sum / sum(slice_ints)) , "]"))
      break
    }
    
    # get the mass of the peak
    slice_mass <- mass_fits[valid_slices[slice_ints[valid_slices] == max(slice_ints[valid_slices])]]
  
    nr_valid_slices <- length(valid_slices)
    valid_slices_array[nr_valid_slices] <- if(is.na(valid_slices_array[nr_valid_slices])) 1 else (valid_slices_array[nr_valid_slices] + 1)
    
    row[(i-1)*2+3] <- slice_mass
    row[(i-1)*2+4] <- slice_ints_sum
  }
    
  if(valid_entry){
    names(row) <- NA
    res_table <- rbind(res_table, row)
    valid_entries <- valid_entries + 1
  }else{
    invalid_entries <- invalid_entries + 1
  }
  
} 

# get the sample names
sample_names <- unlist(lapply(all_datasets, function(d){
  sample_name <- paste(d$sample, d$name, sep=".")
  c(paste(sample_name, "mass", sep="."), paste(sample_name, "int", sep="."))
}))

names(res_table) <- c("protein.ac", "theo.mol.weight", sample_names)

# write thable
write.table(res_table, file = res_path, sep = "\t", row.names = FALSE)
  
