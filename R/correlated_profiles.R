library(httr)
library(RCurl)


rm(list=ls())


# results
res_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/correlated_profiles/pumba_correlated_profiles_220217.txt")

# peaks table
peaks_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/well_behaved/pumba_human_proteins_211118.txt")

# parameters
organism_param <- "human"
db_param <- "UP000005640_canonical_plus_splices_16092020"

# matching peaks to theo mass have a distance less than that (log10)
same_peak_thres <- 0.08
# a strong enough peak is at least this value comparing to strongest peak
peak_int_thres <- 0.1
# peaks that match are pretty precise
matching_peak_thres <- 0.005

#nr_proteins <- 100

# first we look for possible results from the peaks table
peaks_table <- read.table(file=peaks_path, header=TRUE, sep="\t")


# internal library
source("./R/functions/all_functions.R")

all_datasets <- get_all_datasets(organism=organism_param)


# create a list of secondary protein peaks of good matches

potential_peak_list <- list()

for(i in 1:nrow(peaks_table)){
#for(i in 1:100){

  print(paste0(i, " of ", nrow(peaks_table)))
  row <- peaks_table[i,]
  
  for(k in 1: length(all_datasets)){
    #dataset <- paste0(all_datasets[[k]]$sample, ".", all_datasets[[k]]$name)
    dataset <- paste0(all_datasets[[k]]$sample)
    
    protein_ac <- row$protein.ac
    theo_mass <- log10(row$theo.mass)
    
    # keep only first ac
    is_first_ac <- row[paste0(dataset, ".is.first.ac")]
    #if(is.na(is_first_ac) || ! is_first_ac) next
    
    # only keep entries which have more then 1 peak
    nr_peaks <- as.numeric(row[paste0(dataset, ".nr.peaks")])
    if(!is.na(nr_peaks) && nr_peaks < 2) next
    
    mass_list <- strsplit(as.character(row[paste0(dataset, ".peak_masses")]), ";")[[1]]
    mass_list <- unlist(lapply(mass_list, function(x){log10(as.numeric(x))}))
    
    mass_dists <- abs(mass_list - theo_mass)
    
    int_list <- as.numeric(strsplit(as.character(row[paste0(dataset, ".peak_ints")]), ";")[[1]])
    int_perc <- int_list/max(int_list)
    
    # should we remove entries which do not much the theo mass at all ?
    
    # only keep entries with strong enough peaks away from theo mass
    keep_peak_ids <- (int_perc > peak_int_thres & mass_dists > same_peak_thres)
    
    if(is.na(keep_peak_ids) || ! any(keep_peak_ids)) next
    
    new_entries <- data.frame()
    
    for(keep_id in 1:length(keep_peak_ids)){
      if(! keep_peak_ids[keep_id]) next
      new_entry <- data.frame(protein_ac=protein_ac, theo_mass=10^theo_mass, mass=10^mass_list[keep_id], int_perc=int_perc[keep_id])
      new_entries <- rbind(new_entries, new_entry)
    }
    
    if(is.null(potential_peak_list[[dataset]])){
      potential_peak_list[[dataset]] <- new_entries
    }else{
      potential_peak_list[[dataset]] <- rbind(potential_peak_list[[dataset]], new_entries)
    }
  }
}


# check if the searched proteins are inside
potential_peak_list[["K562.12160"]][potential_peak_list[["K562.12160"]]$protein_ac == "Q16777",]



# now we can look in the list for good matching pairs
matches <- list()

for(k in 1: length(all_datasets)){
  dataset <- paste0(all_datasets[[k]]$sample, ".", all_datasets[[k]]$name)
  potential_peaks <- potential_peak_list[[dataset]]
  
  already_taken <- c()
  
   for(i in 1:(nrow(potential_peaks)-1)){
     if(i %in% already_taken) next
     mass <- log10(potential_peaks$mass[i])
     this_matches <- c(potential_peaks$protein_ac[i])
     
     for(j in (i+1):nrow(potential_peaks)){
       if(j %in% already_taken) next
       if(abs(mass - log10(potential_peaks$mass[j])) <= matching_peak_thres){
         this_matches <- c(this_matches, potential_peaks$protein_ac[j])
         already_taken <- c(already_taken, j)
       }
     }
     
     if(length(this_matches) > 1){
       if(is.null(matches[[dataset]])) matches[[dataset]] <- list()
       matches[[dataset]][[length(matches[[dataset]]) + 1]] <- this_matches
     }
     
   }
}


# look for the examples given by Manfredo
# K562.12159 should contain all of them
ex1 <- c("P11274", "P00519")
ex2 <- c("Q16777", "P62979")
ex3 <- c("P63165", "P46060")
ex_dataset <- "K562.12160" #"HeLa.12166" #"K562.12159"

checkEx <- function (ex){
  if(ex[1] %in% a && ex[2] %in% a){
    print(paste0(ex[1], " and ", ex[2], " found"))
    print(a)
    return(1)
  }else{
    return(0)
    #print(paste0(ex1, " or ", ex2, " not found"))
  }
}

for(k in 1: length(all_datasets)){
  ex_dataset <- paste0(all_datasets[[k]]$sample, ".", all_datasets[[k]]$name)
  nr_found <- 0

  for(i in 1:length(matches[[ex_dataset]])){
    a <- matches[[ex_dataset]][[i]]
    nr_found <- nr_found + checkEx(ex1)
    nr_found <- nr_found + checkEx(ex2)
    nr_found <- nr_found + checkEx(ex3)
  }
  
  print(paste0(ex_dataset, " ", nr_found))
}





all_protein_groups <- get_all_protein_groups(all_datasets)

# take all protein_acs
protein_acs <-get_all_protein_acs(all_protein_groups)


# test proteins
protein_acs <- c("Q16777", "P62979", "P63165", "P46060", "P11274",  "P00519")


# loop over all proteins
nr_prot_loop <- if(! exists("nr_proteins")) length(protein_acs) else nr_proteins

