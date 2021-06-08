

# create the cache directory if not already there
seq_data_cache <- paste0(data_cache, "sequences")
if(! dir.exists(seq_data_cache)) dir.create(seq_data_cache)


get_seq <- function(protein_ac, database){
  cache_path <- paste0(seq_data_cache, '/', protein_ac, '.RData')
  
  if(file.exists(cache_path)){
    load(cache_path)
  }else{
    seq_get <- getURL(paste0("localhost:9000/sequence/", protein_ac, "/database/", database))
    seq_object <- RJSONIO::fromJSON(seq_get)
    save(file=cache_path, seq_object)
  }
  seq_object
}


