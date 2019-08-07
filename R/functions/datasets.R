
# get all datasets
get_all_datasets <- function(){
  dataset_get <- GET("localhost:9000/dataset")
  content(dataset_get)
}

# get list of samples
get_samples <- function(all_datasets){
  my_samples <- c()
  for(d in all_datasets){ if(! d$sample %in% my_samples) my_samples <- c(my_samples, d$sample) }
  my_samples
}


# check if the current sample is in the list of selected samples
is_selected_sample <- function(sample, sel_samples){
  length(sel_samples) == 0 || sample %in% sel_samples
}

# find datasets of a given sample
find_datasets <- function(sample, all_datasets){
  my_res <- list()
  i <- 1
  for(d in all_datasets){ if(d$sample == sample){my_res[[i]] <- d; i <- i + 1} }
  my_res
}

# get dataset ids
get_dataset_ids <- function(datasets){
  unlist(lapply(datasets, function(d) d$id ))
}

# get protein groups of a given dataset
get_protein_groups <- function(dataset){
  dataset_pg_path <- paste0(data_root, dataset$massFitResult$proteinGroupsPath)
  read.table(file=dataset_pg_path, sep="\t", header = TRUE)
}

# get the first protein AC from the list
get_protein_acs <- function(protein_groups){
  unlist(lapply(strsplit(as.character(protein_groups$Majority.protein.IDs), ";"), function(x){x[1]})) 
}

# get merged data for a protein from the backend
get_protein_merge <- function(protein_ac, sample){
  # create the folder for the dataset if necessary
  dataset_cache_path <- paste0(data_cache, sample)
  if(! dir.exists(dataset_cache_path)) dir.create(dataset_cache_path)
  
  # load protein from cache or get it from the backend
  protein_cache_path <- paste0(dataset_cache_path, '/', protein_ac, '.RData')
  
  if(file.exists(protein_cache_path)){
    load(protein_cache_path)
  }else{
    protein_merge_get <- GET(paste0("http://localhost:9000/merge-protein/", protein_ac), query=list(dataSetsString=paste(dataset_ids, collapse=",")))
    protein_merge <- content(protein_merge_get)
    save(protein_merge, file=protein_cache_path)
  }
  protein_merge
}

# get masses and intensities of merge
get_masses <- function(protein_merge){ 
  10^(unlist(protein_merge$theoMergedProtein$theoMolWeight))
}

get_ints <- function(protein_merge){
  unlist(protein_merge$theoMergedProtein$intensities)
}



