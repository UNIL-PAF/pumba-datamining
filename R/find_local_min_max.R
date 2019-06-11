library(httr)

data_root <- "/Users/admin/tmp/pumba/"

# get all datasets
dataset_get <- GET("localhost:9000/dataset")
datasets <- content(dataset_get)

# select one dataset
dataset <- datasets[[1]]

for(dataset in datasets){
  # load normalized protein groups file
  dataset_path <- paste0(data_root, dataset$massFitResult$proteinGroupsPath)
  protein_groups <- read.table(file=dataset_path, sep="\t", header = TRUE)
  
  print(dim(protein_groups))
  print(dataset$id)
  
  # select one protein
  protein_ac <- "O00754"
  first_protein_acs <- unlist(lapply(strsplit(as.character(protein_groups$Majority.protein.IDs), ";"), function(x){x[1]}))
  protein <- protein_groups[which(first_protein_acs == protein_ac),]
  
  #protein_ac <- strsplit(as.character(protein$Majority.protein.IDs), ";")[[1]][1]
  
  protein_merge_get <- GET(paste0("http://localhost:9000/merge-protein/", protein_ac), query=list(dataSetsString=dataset$id))
  protein_merge <- content(protein_merge_get)
  
  # for some datasets the protein is not detected
  print(length(protein_merge))
  if(length(protein_merge) == 0) next
  
  # since we load only one dataset at the time there is only 1 list
  protein_merge <- protein_merge[[1]]
  names(protein_merge)
  
  mass <- 10^(unlist(protein_merge$theoMergedProtein$theoMolWeight))
  ints <- unlist(protein_merge$theoMergedProtein$intensities)
  peaks <- data.frame(mass, ints)
  
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
  
  plot(peaks, type="l", xlim=c(10,200))
  #plot(peaks, type="l")
  points(peaks[peaks_idx,], pch=20, col="red")
}


  