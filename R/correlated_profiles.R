library(httr)
library(RCurl)
library(lsa)

rm(list=ls())

# output path
out_path <- "/Users/rmylonas/tmp/datamining_pumba/results/correlated_profiles/"
current_date <- "220224"

# parameters
organism_param <- "human"
db_param <- "UP000005640_canonical_plus_splices_16092020"

# a strong enough peak is at least this value comparing to strongest peak
peak_int_thres <- 0.08
# peaks that match are pretty precise
matching_peak_thres <- 0.008
# similarity of vectors (cosine)
sim_score_thres <- 0.8
# intensity ratio between the 2 proteins
int_ratio_thresh <- 0.1
# the difference between sum of theo_masses and the peak mass should be in this range
theo_mass_range <- c(-5, 70)
max_theo_dist_thresh <- 0.05

#########################
# pathes
#########################
# pumba backend data
data_root <- "/Users/rmylonas/tmp/pumba/"
# cache 
data_cache <- "/Users/rmylonas/tmp/datamining_pumba/"

#########################
# settings for peak detection
#########################
peak_detection_int_thresh <- 1e-08
peak_detection_zero_thresh <- 0.8
peak_detection_min_rel_to_strongest <- 0.01

source("./R/functions/sequences.R")
source("./R/functions/datasets.R")
source("./R/functions/peak_detection.R")

all_datasets <- get_all_datasets(organism=organism_param)

all_protein_groups <- get_all_protein_groups(all_datasets)

# take all protein_acs
protein_acs <-get_all_protein_acs(all_protein_groups)

# test proteins
protein_acs <- c("Q16777", "P62979", "P63165", "P46060", "P11274",  "P00519")

dataset_names <- unlist(lapply(all_datasets, function(x){x$name}))
dataset_samples <- unlist(lapply(all_datasets, function(x){x$sample}))
samples <- unique(dataset_samples)

potential_peak_list <- list()
int_vec_sample <- list()

# loop over all proteins
nr_prot_loop <- if(! exists("nr_proteins")) length(protein_acs) else nr_proteins

for(k in 1:nr_prot_loop){
  print(paste0(k, " of ", nr_prot_loop)) 
  
  protein_ac <- protein_acs[k]
  
  if(length(grep("REV__", protein_ac)) > 0) next
  
  seq_object <- get_seq(protein_ac, db_param)
  theo_mol_weight <- log10(seq_object$molWeight/1000)
  
  sample_protein_merge <- get_protein_merge(protein_ac = protein_ac)
  sample_protein_merge_names <- unlist(lapply(sample_protein_merge$proteinMerges, function(x){ x[['sample']] }))
  protein_name <- sample_protein_merge$mainSequence$proteinName
  
  for(i in 1:length(sample_protein_merge_names)){
    sample_name <- sample_protein_merge_names[i]
    
    sel_theo_merged_prot <- sample_protein_merge$proteinMerges[[i]]$theoMergedProtein
    sel_proteins <- sample_protein_merge$proteinMerges[[i]]$proteins
    sample_peak_idxs <- get_peak_indexes(sel_theo_merged_prot$intensities)
    
    if(is.null(sample_peak_idxs)) next
    
    sample_mass <- sel_theo_merged_prot$theoMolWeights
    sample_peak_masses <- sample_mass[sample_peak_idxs]
    
    # if the theoretical weight has a peak we remove it
    theo_dists <- abs(sample_peak_masses - theo_mol_weight)
    min_theo_dist_idx <- which(theo_dists == min(theo_dists))
    min_theo_dist <- abs(sample_peak_masses[min_theo_dist_idx] - theo_mol_weight)
    
    theo_peak_idx <- if(min_theo_dist <= max_theo_dist_thresh) min_theo_dist_idx else NA

    if(length(sample_peak_idxs) == 0 || (length(sample_peak_idxs) == 1 && !is.na(theo_peak_idx))) next
    
    sample_ints <- sel_theo_merged_prot$intensities
    
    sample_peak_ints <- sample_ints[sample_peak_idxs]
    norm_ints <- sample_peak_ints/max(sample_peak_ints)
    
    #valid_peaks <- (norm_ints >= peak_int_thres)
    
    # only consider highest peak
    peak_int_order <- rank(-1*sample_peak_ints)
    #if(! is.na(theo_peak_idx)) peak_int_order[theo_peak_idx] <- 99
    highest_peak_idx <- c(which(peak_int_order == 1), which(peak_int_order == 2))
    valid_peaks <- rep(FALSE, length(sample_peak_idxs))
    valid_peaks[highest_peak_idx] <- TRUE
    if(! is.na(theo_peak_idx)) valid_peaks[theo_peak_idx] <- FALSE
    
    new_entries <- data.frame()
    int_vec_peak <- list()
    
    # plot(sample_ints)
    # a <- abs(theo_mol_weight-sample_mass)
    # theo_mol_weight_idx <- which(min(a)==a)
    # abline(v=theo_mol_weight_idx, col="green")
    
    for(l in 1:length(valid_peaks)){
      is_valid <- valid_peaks[l]
      if(! is_valid) next
      
      # find the peaks to compute the dot-product
      peak_limits <- get_peak_limits(sample_ints, sample_peak_idxs[l])
      mass_limits <- sample_mass[peak_limits]
      
      #abline(v=sample_peak_idxs[l], col="red", lty=3)
      
      int_vec_rep <- list()
      
      for(m in 1:length(sel_proteins)){
        prot <- sel_proteins[[m]]
        
        p_sample_name <- prot$dataSet$name
        p_mass_fits <- prot$dataSet$massFitResult$massFits
        p_sel_slices <- which(p_mass_fits >= mass_limits[1] & p_mass_fits <= mass_limits[2])
        p_ints <- prot$intensities
        
        int_vec_rep[[p_sample_name]][["ints"]] <- p_ints
        int_vec_rep[[p_sample_name]][["sel_slices"]] <- p_sel_slices
      }
      
      int_vec_peak[[length(int_vec_peak) + 1]] <- int_vec_rep
      
      new_entry <- data.frame(protein_ac=protein_ac, protein_name=protein_name, theo_mass=10^theo_mol_weight, mass=10^sample_peak_masses[l], int=sample_peak_ints[l], int_perc=norm_ints[l])
      new_entries <- rbind(new_entries, new_entry)
    }
    
    int_vec_sample[[sample_name]] <- append(int_vec_sample[[sample_name]], int_vec_peak)
    
    if(is.null(potential_peak_list[[sample_name]])){
      potential_peak_list[[sample_name]][["entries"]] <- new_entries
    }else{
      potential_peak_list[[sample_name]][["entries"]] <- rbind(potential_peak_list[[sample_name]][["entries"]], new_entries)
    }
    
  }
}


similarity_score <- function(sample, idx1, idx2){
  int_vec <- int_vec_sample[[sample]]
  min_val <- 1
  
  for(s_id in 1:length(int_vec[[idx1]])){
    rep_name <-names(int_vec[[idx1]][s_id])
    # return 0 if the samples are not the same
    if(is.na(names(int_vec[[idx2]][s_id])) || names(int_vec[[idx2]][s_id]) != rep_name) return(0)
    
    sel_slices <- intersect(int_vec[[idx1]][[s_id]]$sel_slices, int_vec[[idx2]][[s_id]]$sel_slices)
    a <- int_vec[[idx1]][[s_id]]$ints[sel_slices]
    b <- int_vec[[idx2]][[s_id]]$ints[sel_slices]
    
    if(length(sel_slices) < 3){
      cos_score <- 0
    }else{
      cos_score <- as.numeric(cosine(a, b))
    }
    
    min_val <- min(min_val, cos_score, na.rm = TRUE)
    
  }
  
  return(min_val)
}


# now we can look in the list for good matching pairs
matches <-list()


for(k in 1:length(samples)){
  sample_name <- samples[k]
  potential_peaks <- potential_peak_list[[sample_name]]
  
  this_matches <- data.frame()
  
   for(i in 1:(nrow(potential_peaks$entries)-1)){
     mass <- log10(potential_peaks$entries$mass[i])
     
     for(j in (i+1):nrow(potential_peaks$entries)){
       if(abs(mass - log10(potential_peaks$entries$mass[j])) <= matching_peak_thres){
         sim_score <- similarity_score(sample_name, i, j)
         
         sorted_ints = sort(c(potential_peaks$entries$int[i], potential_peaks$entries$int[j]))
         int_ratio <- sorted_ints[1]/sorted_ints[2]
         theo_mass_diff = potential_peaks$entries$theo_mass[i] + potential_peaks$entries$theo_mass[j] - potential_peaks$entries$mass[i]
         
         if(sim_score >= sim_score_thres && int_ratio >= int_ratio_thresh && theo_mass_diff >= theo_mass_range[1] && theo_mass_diff <= theo_mass_range[2]){
           
           new_row <- data.frame(
             protein_ac_1=potential_peaks$entries$protein_ac[i], 
             protein_ac_2=potential_peaks$entries$protein_ac[j],
             protein_name_1=potential_peaks$entries$protein_name[i],
             protein_name_2=potential_peaks$entries$protein_name[j],
             theo_weight_1=potential_peaks$entries$theo_mass[i],
             theo_weight_2=potential_peaks$entries$theo_mass[j],
             peak_mass_1=potential_peaks$entries$mass[i],
             peak_mass_2=potential_peaks$entries$mass[j],
             similarity_score=sim_score,
             int_ratio=int_ratio,
             theo_mass_diff=theo_mass_diff
             )
           this_matches <- rbind(this_matches, new_row)
         }
       }
     }
   }
  matches[[sample_name]] <- this_matches
}

matches[["K562"]]
matches
