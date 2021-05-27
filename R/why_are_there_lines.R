library(httr)
library(ggplot2)
library(Peptides)

rm(list=ls())

data_root <- "/Users/rmylonas/tmp/pumba/"
data_cache <- "/Users/rmylonas/tmp/datamining_pumba/"

sel_datasets <- c("1559132129057") #NULL #c("1559132129057")
nr_proteins <- 1000

# used to get the charge for each protein
fasta_seq_path <- "/Users/rmylonas/Work/PAF/projects/pumba/data/fasta/20170609_UP000005640_9606.csv"
fasta_seqs <- read.csv(file=fasta_seq_path)

# get all datasets
dataset_get <- GET("localhost:9000/dataset")
datasets <- content(dataset_get)

# store results
results <- list()

for(dataset in datasets){
  # look only at selected datasets
  if(length(sel_datasets) > 0 && ! (dataset$id %in% sel_datasets)) next

    # create the folder for the dataset if necessary
    dataset_name <- paste0(dataset$id, '_', dataset$sample, '_', dataset$name)
    dataset_cache_path <- paste0(data_cache, dataset_name)
    if(! dir.exists(dataset_cache_path)) dir.create(dataset_cache_path)
  
    # load normalized protein groups file
    dataset_pg_path <- paste0(data_root, dataset$massFitResult$proteinGroupsPath)
    protein_groups <- read.table(file=dataset_pg_path, sep="\t", header = TRUE)
    first_protein_acs <- unlist(lapply(strsplit(as.character(protein_groups$Majority.protein.IDs), ";"), function(x){x[1]}))
    
    # results
    closest_peak_dists <- c()
    theo_weights <- c()
    protein_acs <- c()
    closest_peak_masses <- c()
    charges <- c()
    charges_by_length <- c()
    max_ints <- c()
    
    # loop over all proteins
    nr_prot_loop <- if(! exists("nr_proteins")) length(first_protein_acs) else nr_proteins
    for(k in 1:nr_prot_loop){
      print(paste0(dataset$id, ': ', k, " of ", nr_prot_loop))
      
      # select one protein
      protein <- protein_groups[k,]
      prot_ints <- protein[grep("Intensity", colnames(protein))]
      max_int <- which(max(prot_ints) == prot_ints)
      protein_ac <- first_protein_acs[k]
      
      # load protein from cache or get it from the backend
      protein_cache_path <- paste0(dataset_cache_path, '/', protein_ac, '.RData')
      if(exists("protein_merge")) rm(protein_merge) 
      
      if(file.exists(protein_cache_path)){
        load(protein_cache_path)
      }else{
        protein_merge_get <- GET(paste0("http://localhost:9000/merge-protein/", protein_ac), query=list(dataSetsString=dataset$id))
        protein_merge <- content(protein_merge_get)
        save(protein_merge, file=protein_cache_path)
      }
      
      # for some datasets the protein is not detected
      if(length(protein_merge) == 0) next
      
      # since we load only one dataset at the time there is only 1 list
      protein_merge <- protein_merge[[1]]
      names(protein_merge)
      
      mass <- 10^(unlist(protein_merge$theoMergedProtein$theoMolWeight))
      ints <- unlist(protein_merge$theoMergedProtein$intensities)  
      
      int_thresh <- 1e-07
      zero_thresh <- 1e-8
      
      i <- 2
      # test if we're in a increasing part of the curve
      is_increasing <- FALSE
      # assure that we passed close to 0 before looking the next curve
      passed_by_zero <- TRUE
      # if neighboring peaks don't pass at 0 in between, we take the highest one.
      last_int <- 0
      peaks_idx <- c()
      
      while(i <= length(ints)){
        if(ints[i] > ints[i-1]){
          is_increasing <- TRUE
          if(! passed_by_zero && ints[i-1] <= zero_thresh) passed_by_zero <- TRUE
        }else{
          if(is_increasing){
            if(ints[i-1] > int_thresh){
              if(passed_by_zero){
                peaks_idx <- c(peaks_idx, i-1)
                passed_by_zero <- FALSE
                last_int <- ints[i-1]
              }else if(ints[i-1] > last_int){
                last_int <- ints[i-1]
                peaks_idx[length(peaks_idx)] <- i-1
              }
            }
            is_increasing <- FALSE
          }
        }
        i <- i+1
      }
      
      # just skip with a warning if no peak was detected
      if(length(peaks_idx) < 1){
        warning(paste0("Could not find any peak in [", protein_ac, "]."))
        next
      }
      
      # filter out peaks below certain threshold relative to strongest peak
      min_rel_to_strongest <- 0.1
      peaks_ints <- ints[peaks_idx]
      min_thresh <- max(peaks_ints) * min_rel_to_strongest
      peaks_ints_ok <- peaks_ints >= min_thresh
      peaks_idx <- peaks_idx[peaks_ints_ok]
      
      peaks_masses <- mass[peaks_idx]
      theo_weight <- as.numeric(protein$Mol..weight..kDa.)
      peak_dists <- abs(peaks_masses - theo_weight)
      
      closest_peak <- which(peak_dists == min(peak_dists))
      closest_peak_perc_diff <- (peaks_masses[closest_peak] / theo_weight) - 1
      peak_dists_perc <- (peaks_masses / theo_weight) - 1
      # peak_dists_log <- log(peaks_masses) - log(theo_weight)
      # closest_peak_log_diff <- log(peaks_masses[closest_peak]) - log(theo_weight)
      
      if(sum(fasta_seqs$ac == protein_ac) == 1){
        seq <- as.character(fasta_seqs$seq[fasta_seqs$ac == protein_ac])
        
        theo_weights <- c(theo_weights, theo_weight)
        closest_peak_dists <- c(closest_peak_dists, closest_peak_perc_diff)
        protein_acs <- c(protein_acs, protein_ac)
        closest_peak_masses <- c(closest_peak_masses, peaks_masses[closest_peak])
        this_charge <- charge(seq, pH=7, pKscale="EMBOSS")
        charges <- c(charges, this_charge)
        charges_by_length <- c(charges_by_length, this_charge / nchar(seq))
        max_ints <- c(max_ints, max_int)
      }
      
      # plot(mass, ints, type="l", main=protein_ac)
      # abline(v=theo_weight, col="blue")
      # points(mass[peaks_idx], ints[peaks_idx], col="red")
      # text(mass[peaks_idx], ints[peaks_idx], labels=(round(peak_dists_log, digits = 2)), col="red", pos=4)
    }
    
    peak_dists <- data.frame(theo_weights, closest_peak_dists, protein_acs, closest_peak_masses, charges, charges_by_length, max_ints)
    peak_dists$color <- "neutral"
    peak_dists$color[peak_dists$max_ints == 44] <- "pos"
    #peak_dists$color[peak_dists$charges < -5] <- "neg"
    
    # peak_dists$color <- "neutral"
    # peak_dists$color[peak_dists$charges_by_length > 0.02] <- "pos"
    # peak_dists$color[peak_dists$charges_by_length < -0.02] <- "neg"
    
    p <- ggplot(data=peak_dists, aes(theo_weights, closest_peak_dists, colour=color)) 
    p <- p + geom_point(alpha=0.2)
    #p <- p + geom_text(data=peak_dists[peak_dists$theo_weights > 200, ], aes(theo_weights, closest_peak_dists, label=protein_acs), srt=90, col="red", alpha=0.3)
    p <- p + xlim(0, 100)
    p <- p + ylim(-2, 2)
    p <- p + scale_colour_manual(values=c("pos" = "green", "neg" = "red", "neutral" = "grey"))
    p <- p + geom_hline(yintercept = 0)
    print(p)
    
    # store the data in the results list
    results[[dataset_name]] <- peak_dists
}
  
# perc_diff = (obs / theo) - 1

  # protein_ac <- "O00754"
  # protein <- protein_groups[which(first_protein_acs == protein_ac),]
  # protein_ac <- strsplit(as.character(protein$Majority.protein.IDs), ";")[[1]][1]
  
 
#   peaks <- data.frame(mass, ints)
#   
# 
#   plot(peaks, type="l", xlim=c(10,200))
#   #plot(peaks, type="l")
#   points(peaks[peaks_idx,], pch=20, col="red")
#   
#   peaks_masses <- mass[peaks_idx]
# #}
# 
# 
# 



# 
# plot(theo_weights, closest_peak_dists, xlim=c(1,500))
# 
# gt_200 <- theo_weights > 200
# 
# 
# text(theo_weights[gt_200], closest_peak_dists[gt_200], labels = protein_acs[gt_200], col="red", srt=90)
# 
# 
# 
# 
# 
# ac_row <- c("Q9Y4A5", "Q03164", "Q709C8", "Q76CQ2", "P46939", "Q14789", "P49454")
# a <- which(first_protein_acs %in% ac_row)
# protein_groups[a, ]
# 
# # find which column provides max value
# int_cols <- grep("Intensity", colnames(protein_groups))
# max_int_slices <- unlist(apply(protein_groups[, int_cols], 1, function(x){which(x == max(x))}))
# map_prot_ids <- which(first_protein_acs %in% protein_acs)
# max_int_slices_2 <- max_int_slices[map_prot_ids]








