library(httr)
library(RCurl)
library(RJSONIO)
library(ggplot2)
library(Peptides)

rm(list=ls())

# results
res_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/well_behaved/pumba_human_proteins_210817.txt")

# parameters
organism_param <- "human"
db_param <- "UP000005640_canonical_plus_splices_16092020"

# internal library
source("./R/functions/all_functions.R")

all_datasets <- get_all_datasets(organism=organism_param)

#nr_proteins <- 10
#plot_proteins <- c("A0A096LP01", "P0DPI2", "A0A0B4J2D5")

all_protein_groups <- get_all_protein_groups(all_datasets)

# just take the protein_acs from the first dataset
# protein_acs <- get_protein_acs(all_protein_groups[[1]])

# take all protein_acs
protein_acs <-get_all_protein_acs(all_protein_groups)

# loop over all proteins
nr_prot_loop <- if(! exists("nr_proteins")) length(protein_acs) else nr_proteins

dataset_ids <- get_dataset_ids(all_datasets)

# prepare the results table
# fields:
#
dataset_names <- unlist(lapply(all_datasets, function(x){x$name}))
dataset_samples <- unlist(lapply(all_datasets, function(x){x$sample}))
samples <- unique(dataset_samples)

column_field_names <- c("is.first.ac", "nr.peaks",
                        "nr.peptides", "peak_masses", "peak_ints", "secondary_peak_masses", "highest.is.closest", 
                        "highest.peak.int", "highest.peak.int.perc", "highest.peak.mass",
                        "highest.peak.start", "highest.peak.end", "closest.peak.int",
                        "closest.peak.int.perc", "closest.peak.mass", "closest.peak.start",
                        "closest.peak.end")

sample_field_names <- c("highest.peak.int", "highest.peak.int.perc", "highest.peak.mass")

column_names <- unlist(lapply(samples, function(sample_name){
  sub_column_names <- paste0(sample_name, ".", sample_field_names)
  dataset_subset <- dataset_samples == sample_name
  column_sample_names <- paste0(dataset_samples[dataset_subset], ".", dataset_names[dataset_subset])    
  c(sub_column_names, unlist(lapply(column_sample_names, function(x){paste0(x, ".", column_field_names)})))
}))

column_field_size <- length(column_field_names)
sample_field_size <- length(sample_field_names)

column_names <- c("protein.ac", "theo.mass", "pI", "hydrophobicity", column_names)
res_table = data.frame(matrix(ncol=length(column_names),nrow=0, dimnames=list(NULL, column_names)))

for(k in 1:nr_prot_loop){
  
  print(paste0(k, " of ", nr_prot_loop))
  protein_ac <- protein_acs[k]
  
  # skip proteins containing the term "REV"
  if(grepl("REV", protein_ac, fixed=TRUE)) next
  
  seq_object <- get_seq(protein_ac, db_param)
  theo_mol_weight <- seq_object$molWeight/1000
  
  # sequence properties
  my_pI <- pI(seq_object$sequence, pKscale = "EMBOSS")
  my_hydrophopicity <- hydrophobicity(seq_object$sequence, scale = "HoppWoods")
  
  one_row = data.frame(matrix(ncol=length(column_names),nrow=1, dimnames=list(NULL, column_names)))
  one_row[1] <- protein_ac
  one_row[2] <- theo_mol_weight
  one_row[3] <- my_pI
  one_row[4] <- my_hydrophopicity
  
  sample_protein_merge <- get_protein_merge(protein_ac = protein_ac)
  sample_protein_merge_names <- unlist(lapply(sample_protein_merge$proteinMerges, function(x){ x[['sample']] }))
  
  for(j in 1:length(samples)){
    
    sub_dataset_ids <- dataset_ids[dataset_samples == samples[j]]
    column_j <- if(exists("column_i")) {column_i + column_field_size} else { sample_field_size * (j-1) + 5 }
    
    sel_protein_merge_id <- which(sample_protein_merge_names == samples[j])
    
    # some samples don't have no results
    if(length(sel_protein_merge_id) == 0){
      next
    }
    
    sel_theo_merged_prot <- sample_protein_merge$proteinMerges[[sel_protein_merge_id]]$theoMergedProtein
    
    # get the indexes of intensities corresponding to a peak (see "config.R" for thresholds)
    sample_peak_idxs <- get_peak_indexes(sel_theo_merged_prot$intensities)
    
    # only fill the columns if peaks where found (otherwise it stays NA's)
    if(length(sample_peak_idxs) >= 1){
      sample_ints <- sel_theo_merged_prot$intensities
      sample_mass <- sel_theo_merged_prot$theoMolWeights
      sample_peak_ints <- sample_ints[sample_peak_idxs]
      sample_highest_peak_idx <- sample_peak_idxs[which(sample_peak_ints == max(sample_peak_ints, na.rm=TRUE))]
      sample_highest_int <- sample_ints[sample_highest_peak_idx]
      sample_highest_mass <- sample_mass[sample_highest_peak_idx]
      
      # compute peak percentage
      sample_highest_peak_limits <- get_peak_limits(sample_ints, sample_highest_peak_idx)
      sample_highest_valid_slices <- which(sample_mass >= sample_mass[sample_highest_peak_limits[1]] & sample_mass <= sample_mass[sample_highest_peak_limits[2]])
      sample_highest_slice_ints_sum <- sum(sample_ints[sample_highest_valid_slices])
      sample_highest_peak_int_perc <- sample_highest_slice_ints_sum / sum(sample_ints)
      
      one_row[column_j] <- sample_highest_int
      one_row[column_j + 1] <- sample_highest_peak_int_perc
      one_row[column_j + 2] <- 10^sample_highest_mass
    }
  
    for(i in 1:length(sub_dataset_ids)){
  
      # column index in the row
      column_i <- column_field_size * (i-1) + column_j + sample_field_size
      
      dataset_id <- dataset_ids[i]
      sample_name <- all_datasets[[i]]$sample
      repl_name <- all_datasets[[i]]$name
      
      # get merged data from backend or cache
      protein_merges <- get_single_protein_merge(protein_ac, dataset_id)
      
      # for some datasets the protein is not detected
      if(length(protein_merges$proteinMerges) == 0){
        next
      }
      
      # since we load only one dataset at the time there is only 1 list
      protein_merge <- protein_merges$proteinMerges[[1]]
      
      if(length(protein_merge$proteins) > 1){
        print("multiple proteins!!!")
        print(protein_ac)
        print(length(protein_merge$proteins))
        cat(file="/Users/rmylonas/tmp/multi_proteins.txt", append=TRUE, paste0(protein_ac, ", ", repl_name, "\n"))
      }
      
      mass <- get_masses(protein_merge)
      ints <- get_ints(protein_merge)
      
      # get the indexes of intensities corresponding to a peak (see "config.R" for thresholds)
      peak_idxs <- get_peak_indexes(ints)
      
      # only fill the columns if peaks where found (otherwise it stays NA's)
      if(length(peak_idxs) >= 1){
        if(grep(protein_ac, protein_merge$proteins[[1]]$proteinIDs) > 1){
          is_first_ac <- FALSE
        }else{
          is_first_ac <- TRUE
        }
        one_row[column_i] <- is_first_ac
        
        nr_of_peaks <- length(peak_idxs)
        one_row[column_i + 1] <- nr_of_peaks
        
        nr_peptides <- length(protein_merge$proteins[[1]]$peptides)
        one_row[column_i + 2] <- nr_peptides
        
        peaks_masses <- mass[peak_idxs]
        peaks_ints <- ints[peak_idxs]
        
        one_row[column_i + 3] <- paste(round(peaks_masses, digits=2), collapse = ";")
        one_row[column_i + 4] <- paste(peaks_ints, collapse = ";")
        
        highest_peak <- which(peaks_ints == max(peaks_ints, na.rm=TRUE))
        highest_peak_idx <- peak_idxs[highest_peak]
        highest_peak_limits <- get_peak_limits(ints, highest_peak_idx)
        highest_peak_limits_masses <- mass[highest_peak_limits]
        highest_peak_int <- peaks_ints[highest_peak]
        highest_peak_mass <- peaks_masses[highest_peak]
        highest_peak_dist <- abs(theo_mol_weight - highest_peak_mass)
        
        secondary_peak_masses <- peaks_masses[-highest_peak]
        one_row[column_i + 5] <- paste(round(secondary_peak_masses, digits=2), collapse = ";")
        
        # find the correct limits in the slices
        log_mass_fits <- all_datasets[[i]]$massFitResult$massFits  
        mass_fits <- 10^(unlist(log_mass_fits))
        
        slice_ints <- unlist(protein_merge$proteins[[1]]$intensities)
        highest_valid_slices <- which(mass_fits >= highest_peak_limits_masses[1] & mass_fits <= highest_peak_limits_masses[2])
        highest_slice_ints_sum <- sum(slice_ints[highest_valid_slices])
        highest_peak_int_perc <- highest_slice_ints_sum / sum(slice_ints)
        
        one_row[column_i + 7] <- highest_peak_int
        one_row[column_i + 8] <- highest_peak_int_perc
        one_row[column_i + 9] <- highest_peak_mass
        one_row[column_i + 10] <- highest_peak_limits_masses[1]
        one_row[column_i + 11] <- highest_peak_limits_masses[2]
        
        peak_theo_dists <- abs(theo_mol_weight - peaks_masses)
        closest_peak <- which(peak_theo_dists == min(peak_theo_dists, na.rm=TRUE))
        closest_peak_idx <- peak_idxs[closest_peak]
        closest_peak_limits <- get_peak_limits(ints, closest_peak_idx)
        closest_peak_limits_masses <- mass[closest_peak_limits]
        closest_peak_int <- peaks_ints[closest_peak]
        closest_peak_mass <- peaks_masses[closest_peak]
        closest_peak_dist <- abs(theo_mol_weight - closest_peak_mass)
        
        highest_is_closest <- highest_peak == closest_peak
        one_row[column_i + 6] <- highest_is_closest
        
        one_row[column_i + 12] <- closest_peak_int
        
        closest_valid_slices <- which(mass_fits >= closest_peak_limits_masses[1] & mass_fits <= closest_peak_limits_masses[2])
        closest_slice_ints_sum <- sum(slice_ints[closest_valid_slices])
        closest_peak_int_perc <- closest_slice_ints_sum / sum(slice_ints)
        one_row[column_i + 13] <- closest_peak_int_perc
        one_row[column_i + 14] <- closest_peak_mass
        one_row[column_i + 15] <- closest_peak_limits_masses[1]
        one_row[column_i + 16] <- closest_peak_limits_masses[2]
      }
      
    }
    
  }
  
  rm("column_i")
  res_table <- rbind(res_table, one_row)
  
  # if(protein_ac %in% plot_proteins){
  #   ## plot the merge curve, the peaks and the distances
  #   peak_dists_perc <- (peaks_masses / theo_mol_weight) - 1
  #   plot(mass, ints, type="l", main=protein_ac)
  #   abline(v=theo_mol_weight, col="blue")
  #   points(mass[peak_idxs], ints[peak_idxs], col="red")
  #   #text(mass[peak_idxs], ints[peak_idxs], labels=(round(peak_dists_perc, digits = 2)), col="red", pos=4)
  # }
  
} 

# write thable
write.table(res_table, file = res_path, sep = "\t", row.names = FALSE)
  
