library(httr)
library(ggplot2)
library(Peptides)

rm(list=ls())

# results
res_path <- ("/Users/admin/tmp/datamining_pumba/results/")
res_name <- "test_"


# internal library
source("./R/functions/all_functions.R")

all_datasets <- get_all_datasets()
samples <- get_samples(all_datasets)

sel_samples <- c("HCT") # NULL
#nr_proteins <- 100

# store results
results <- list()

for(sample in samples){
    # look only at selected datasets
    if(! is_selected_sample(sample, sel_samples)) next
  
    # load normalized protein groups file of one dataset
    datasets <- find_datasets(sample, all_datasets)
    
    # get the dataset ids for this sample
    dataset_ids <- get_dataset_ids(datasets)
    
    protein_groups <- get_protein_groups(datasets[[1]])
    first_protein_acs <- get_protein_acs(protein_groups)
    
    # results
    closest_peak_dists <- c()
    theo_weights <- c()
    protein_acs <- c()
    closest_peak_masses <- c()
    charges <- c()
    charges_by_length <- c()
    pIs <- c()
    locations <- c()
    
    # loop over all proteins
    nr_prot_loop <- if(! exists("nr_proteins")) length(first_protein_acs) else nr_proteins
    for(k in 1:nr_prot_loop){
      print(paste0(sample, ': ', k, " of ", nr_prot_loop))
      
      # select one protein
      protein <- protein_groups[k,]
      protein_ac <- first_protein_acs[k]
      
      # get merged data from backend or cache
      protein_merges <- get_protein_merge(protein_ac, sample)
      
      # for some datasets the protein is not detected
      if(length(protein_merges) == 0) next
      
      # since we load only one dataset at the time there is only 1 list
      protein_merge <- protein_merges[[1]]
      
      mass <- get_masses(protein_merge)
      ints <- get_ints(protein_merge)
      
      # get the indexes of intensities corresponding to a peak (see "config.R" for thresholds)
      peak_idxs <- get_peak_indexes(ints)
      
      # just skip with a warning if no peak was detected
      if(length(peak_idxs) < 1){
        warning(paste0("Could not find any peak in [", protein_ac, "]."))
        next
      }
      
      peaks_masses <- mass[peak_idxs]
      theo_weight <- as.numeric(protein$Mol..weight..kDa.)
      peak_dists <- abs(peaks_masses - theo_weight)
      
      closest_peak <- which(peak_dists == min(peak_dists))
      closest_peak_perc_diff <- (peaks_masses[closest_peak] / theo_weight) - 1
      
      # peak_dists_log <- log(peaks_masses) - log(theo_weight)
      # closest_peak_log_diff <- log(peaks_masses[closest_peak]) - log(theo_weight)
      
      if(seq_exists(protein_ac)){
        seq <- get_seq(protein_ac)
        theo_weights <- c(theo_weights, theo_weight)
        protein_acs <- c(protein_acs, protein_ac)

        # distance info        
        closest_peak_dists <- c(closest_peak_dists, closest_peak_perc_diff)
        closest_peak_masses <- c(closest_peak_masses, peaks_masses[closest_peak])

        # info from AA sequence
        this_charge <- charge(seq, pH=7, pKscale="EMBOSS")
        this_pI <- pI(seq, pKscale = "EMBOSS")
        pIs <- c(pIs, this_pI)
        charges <- c(charges, this_charge)
        charges_by_length <- c(charges_by_length, this_charge / nchar(seq))        
        
        # info from uniprot
        uniprot_xml <- get_uniprot_xml(protein_ac)
        location <- NA
        if(! is.null(uniprot_xml)){
          location <- paste(get_locations(uniprot_xml), collapse=",")
        }
        locations <- c(locations, location)
      }
      
      ## plot the merge curve, the peaks and the distances
      # peak_dists_perc <- (peaks_masses / theo_weight) - 1
      # plot(mass, ints, type="l", main=protein_ac)
      # abline(v=theo_weight, col="blue")
      # points(mass[peaks_idx], ints[peaks_idx], col="red")
      # text(mass[peaks_idx], ints[peaks_idx], labels=(round(peak_dists_perc, digits = 2)), col="red", pos=4)
    }
    
    peak_dists <- data.frame(
      theo_weights, 
      closest_peak_dists, 
      protein_acs, 
      closest_peak_masses, 
      charges, 
      charges_by_length, 
      pIs,
      locations)
    
    # store the data in the results list
    results[[sample]] <- peak_dists
}


# store results with timestamp
res_file <- paste0(res_path, res_name, as.numeric(Sys.time()), ".RData")
save(results, file=res_file)
print(paste0("saved results in [", res_file, "]"))



