library(xml2)

# create the cache directory if not already there
uniprot_data_cache <- paste0(data_cache, "uniprot")
if(! dir.exists(uniprot_data_cache)) dir.create(uniprot_data_cache)

# get the uniprot xml for a given protein_ac either from cache or the web
get_uniprot_xml <- function(protein_ac){
  cache_path <- paste0(uniprot_data_cache, '/', protein_ac, '.xml')
  xml_data <- NULL
  
  if(file.exists(cache_path)){
    xml_data <- xml_ns_strip(read_xml(cache_path))
  }else{
    url_path <- paste0(uniprot_url, protein_ac, uniprot_url_suffix)
    raw_xml_data <- tryCatch(read_xml(url_path), 
                             error=function(cond){
                               message(paste0("Could not load [", protein_ac, "] from [",url_path, "]"))
                               return(NULL)
                             })
    if(! is.null(raw_xml_data)){
      xml_data <- xml_ns_strip(raw_xml_data)
      write_xml(xml_data, file=cache_path) 
    }
  }
  xml_data
}


# get protein locations (e.g. "Cell membrane")
get_locations <- function(uniprot_xml){
  locations <- xml_find_all(uniprot_xml, ".//location")
  locations <- locations[! is.na(xml_attr(locations, "evidence"))]
  unique(xml_text(locations))
}
