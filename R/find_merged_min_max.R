library(httr)
library(ggplot2)
library(Peptides)

rm(list=ls())

# results
res_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/")
res_name <- "homodimers_"

# internal library
source("./R/functions/all_functions.R")

all_datasets <- get_all_datasets()
samples <- get_samples(all_datasets)

sel_samples <- c("HCT") # NULL
#nr_proteins <- 100
plot_proteins <- c() #"H7C417", "O60307", "E9PLN8", "E7EVH7", "A2RTX5", "A2RRP1", "A0A0A6YYH1")

# store results
results <- list()

for(sample in samples){
    # look only at selected datasets
    if(! is_selected_sample(sample, sel_samples)) next
  
    # load normalized protein groups file of one dataset
    datasets <- find_datasets(sample, all_datasets)
    
    # get the dataset ids for this sample
    dataset_ids <- get_dataset_ids(datasets)
    
    all_protein_groups <- get_all_protein_groups(datasets)
    #orig_protein_groups <- get_orig_protein_groups(datasets[[1]])
    all_orig_protein_groups <- get_all_orig_protein_groups(datasets)
    first_protein_acs <- get_protein_acs(all_protein_groups[[1]])
    
    # results
    closest_peak_dists <- c()
    theo_weights <- c()
    protein_acs <- c()
    closest_peak_masses <- c()
    highest_peak_masses <- c()
    charges <- c()
    charges_by_length <- c()
    pIs <- c()
    locations <- c()
    glycosylations <- c()
    signal_peps <- c()
    ptms <- c()
    prot_intensities <- c()
    nr_peaks <- c()
    perc_slices <- c()
    perc_dists <- c()
    peak_intensities <- c()
    peak_masses <- c()
    scop_classes <- c()
    hydrophobicities <- c()
    homodimers <- c()
    
    # loop over all proteins
    nr_prot_loop <- if(! exists("nr_proteins")) length(first_protein_acs) else nr_proteins
    for(k in 1:nr_prot_loop){
      print(paste0(sample, ': ', k, " of ", nr_prot_loop))
      
      # select one protein
      protein <- all_protein_groups[[1]][k,]
      protein_ac <- first_protein_acs[k]
      #orig_protein <- orig_protein_groups[k,]
      
      # look for percentage of slices in which there was a signal
      perc_slices_prot <- get_perc_slices(all_protein_groups, protein_ac, 0.1)
      
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
      peaks_ints <- ints[peak_idxs]
      theo_weight <- as.numeric(protein$Mol..weight..kDa.)
      peak_dists <- abs(peaks_masses - theo_weight)
      peak_perc_dists <- (peaks_masses / theo_weight) - 1
      
      closest_peak <- which(peak_dists == min(peak_dists))
      closest_peak_perc_diff <- (peaks_masses[closest_peak] / theo_weight) - 1
      highest_peak <- which(peaks_ints == max(peaks_ints))
      
      # peak_dists_log <- log(peaks_masses) - log(theo_weight)
      # closest_peak_log_diff <- log(peaks_masses[closest_peak]) - log(theo_weight)
      
      if(seq_exists(protein_ac)){
        seq <- get_seq(protein_ac)
        theo_weights <- c(theo_weights, theo_weight)
        protein_acs <- c(protein_acs, protein_ac)

        # distance info        
        closest_peak_dists <- c(closest_peak_dists, closest_peak_perc_diff)
        closest_peak_masses <- c(closest_peak_masses, peaks_masses[closest_peak])
        highest_peak_masses <- c(highest_peak_masses, peaks_masses[highest_peak])
        
        # number of slices or peaks per protein
        nr_peaks <- c(nr_peaks, length(peaks_masses))
        perc_slices <- c(perc_slices, perc_slices_prot)

        # info from AA sequence
        this_charge <- charge(seq, pH=7, pKscale="EMBOSS")
        this_pI <- pI(seq, pKscale = "EMBOSS")
        this_hydrophopicity <- hydrophobicity(seq, scale = "HoppWoods")
        pIs <- c(pIs, this_pI)
        charges <- c(charges, this_charge)
        charges_by_length <- c(charges_by_length, this_charge / nchar(seq))
        hydrophobicities <- c(hydrophobicities, this_hydrophopicity)
        
        # info from MaxQuan
        prot_intensities <- c(prot_intensities, get_median_protein_intensities(all_orig_protein_groups, protein_ac))
        
        # info from uniprot
        uniprot_xml <- get_uniprot_xml(protein_ac)
        location <- NA
        glycosylation <- NA
        signal_pep <- NA
        ptm <- NA
        is_homodimer <- NA
        if(! is.null(uniprot_xml)){
          location <- paste(get_locations(uniprot_xml), collapse=",")
          glycosylation <- paste(get_glycosylations(uniprot_xml), collapse = ",")
          signal_pep <- get_signal_pep(uniprot_xml)
          ptm <- paste(get_ptms(uniprot_xml), get_crosslink(uniprot_xml), collapse = ",")
          is_homodimer <- get_homodimer(uniprot_xml)
        }
        locations <- c(locations, location)
        glycosylations <- c(glycosylations, glycosylation)
        signal_peps <- c(signal_peps, signal_pep)
        ptms <- c(ptms, ptm)
        peak_intensities <- c(peak_intensities, paste(peaks_ints, collapse = ","))
        peak_masses <- c(peak_masses, paste(peaks_masses, collapse = ","))
        perc_dists <- c(perc_dists, paste(peak_perc_dists, collapse=","))
        scop_classes <- c(scop_classes, get_scop_classes(protein_ac))
        homodimers <- c(homodimers, is_homodimer)
      }
      
      if(protein_ac %in% plot_proteins){
        ## plot the merge curve, the peaks and the distances
        peak_dists_perc <- (peaks_masses / theo_weight) - 1
        plot(mass, ints, type="l", main=protein_ac)
        abline(v=theo_weight, col="blue")
        points(mass[peak_idxs], ints[peak_idxs], col="red")
        #text(mass[peak_idxs], ints[peak_idxs], labels=(round(peak_dists_perc, digits = 2)), col="red", pos=4)
      }
 
    }
    
    sample_res <- data.frame(
      theo_weights, 
      closest_peak_dists, 
      protein_acs, 
      closest_peak_masses,
      highest_peak_masses,
      charges, 
      charges_by_length, 
      pIs,
      locations, 
      glycosylations,
      signal_peps,
      ptms,
      prot_intensities,
      nr_peaks,
      perc_slices,
      peak_intensities,
      peak_masses,
      perc_dists, 
      scop_classes,
      hydrophobicities,
      homodimers
      )
    
    # store the data in the results list
    results[[sample]] <- sample_res
}

# store results with timestamp
res_file <- paste0(res_path, res_name, as.numeric(Sys.time()), ".RData")
save(results, file=res_file)
print(paste0("saved results in [", res_file, "]"))

