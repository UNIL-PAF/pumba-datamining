library(httr)
library(RCurl)
library(lsa)


rm(list=ls())

# output path
out_path <- "/Users/rmylonas/tmp/datamining_pumba/results/correlated_profiles/"
current_date <- "220223"

# parameters
organism_param <- "human"
db_param <- "UP000005640_canonical_plus_splices_16092020"

# matching peaks to theo mass have a distance less than that (log10)
same_peak_thres <- 0.08
# a strong enough peak is at least this value comparing to strongest peak
peak_int_thres <- 0.08
# peaks that match are pretty precise
matching_peak_thres <- 0.01
# similarity of vectors (cosine)
sim_score_thres <- 0.8

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
peak_detection_int_thresh <- 1e-06
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
#protein_acs <- c("Q16777", "P62979", "P63165", "P46060", "P11274",  "P00519")

# loop over all proteins
nr_prot_loop <- if(! exists("nr_proteins")) length(protein_acs) else nr_proteins

dataset_names <- unlist(lapply(all_datasets, function(x){x$name}))
dataset_samples <- unlist(lapply(all_datasets, function(x){x$sample}))
samples <- unique(dataset_samples)

potential_peak_list <- list()
int_vec_sample <- list()

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
    
    if(length(sample_peak_idxs) >= 2){
      sample_ints <- sel_theo_merged_prot$intensities
      sample_mass <- sel_theo_merged_prot$theoMolWeights
      
      sample_peak_ints <- sample_ints[sample_peak_idxs]
      norm_ints <- sample_peak_ints/max(sample_peak_ints)
      sample_peak_masses <- sample_mass[sample_peak_idxs]
      
      valid_peaks <- (abs(sample_peak_masses - theo_mol_weight) >= same_peak_thres & norm_ints >= peak_int_thres)
      
      new_entries <- data.frame()
      int_vec_peak <- list()
      
      for(l in 1:length(valid_peaks)){
        is_valid <- valid_peaks[l]
        if(! is_valid) next
        
        # find the peaks to compute the dot-product
        peak_limits <- get_peak_limits(sample_ints, sample_peak_idxs[l])
        mass_limits <- sample_mass[peak_limits]
        
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
        
        new_entry <- data.frame(protein_ac=protein_ac, protein_name=protein_name, theo_mass=10^theo_mol_weight, mass=10^sample_peak_masses[l], int_perc=norm_ints[l])
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
matches <- list()

for(k in 1:length(samples)){
  sample_name <- samples[k]
  potential_peaks <- potential_peak_list[[sample_name]]
  
  already_taken <- c()
  
   for(i in 1:(nrow(potential_peaks$entries)-1)){
     if(i %in% already_taken) next
     mass <- log10(potential_peaks$entries$mass[i])
     this_matches <- list()
     this_matches[["protein_ac"]] <- c(potential_peaks$entries$protein_ac[i])
     this_matches[["protein_name"]] <- c(potential_peaks$entries$protein_name[i])
     this_matches[["theo_mass"]] <- c(potential_peaks$entries$theo_mass[i])
     
     for(j in (i+1):nrow(potential_peaks$entries)){
       if(j %in% already_taken) next
       if(abs(mass - log10(potential_peaks$entries$mass[j])) <= matching_peak_thres){
         already_taken <- c(already_taken, j)
         sim_score <- similarity_score(sample_name, i, j)

         if(sim_score >= sim_score_thres){
           this_matches[["protein_ac"]] <- c( this_matches[["protein_ac"]], potential_peaks$entries$protein_ac[j])
           this_matches[["protein_name"]] <- c( this_matches[["protein_name"]], potential_peaks$entries$protein_name[j])
           this_matches[["theo_mass"]] <- c( this_matches[["theo_mass"]], potential_peaks$entries$theo_mass[j])
         }
       }
     }
     
     if(length(this_matches$protein_ac) > 1){
       this_matches[["mass"]] <- 10^mass
       if(is.null(matches[[sample_name]])) matches[[sample_name]] <- list()
       matches[[sample_name]][[length(matches[[sample_name]]) + 1]] <- this_matches
     }
     
   }
}

matches


# we create a file for each sample

for(k in 1:length(samples)){
  sample_name <- samples[k]
  res_file <- paste0(out_path, "corr_profiles_", sample_name, "_", current_date, ".txt")
  one_match <- matches[[sample_name]]
  
  cat("", file=res_file)
  
  for(l in 1:length(one_match)){
    m <- one_match[[l]]
    cat(l, ":\t", m$mass, "\n",  file=res_file, append=TRUE)
    cat(paste0(m$protein_ac, " (", m$protein_name, ") [", m$theo_mass, "]", collapse="; "), "\n", file=res_file, append=TRUE)
    cat("\n", file=res_file, append = TRUE)
  }
}



