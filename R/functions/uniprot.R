library(xml2)

# valid ubiquitins
valid_ubiquitin <- c("SUMO1", "SUMO2" , "SUMO3", "ubiquitin", "NEDD8", "ISG15", "ATG5", "KEAP1", "ATG12", "UFM1", "SUMO")

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
    print(paste0("Download ", protein_ac))
    url_path <- paste0(uniprot_url, protein_ac, uniprot_url_suffix)
    raw_xml_data <- tryCatch(read_xml(url_path), 
                             error=function(cond){
                               message(paste0("Could not load [", protein_ac, "] from [",url_path, "]"))
                               return(NULL)
                             })
    if(! is.null(raw_xml_data)){
      xml_data <- xml_ns_strip(raw_xml_data)
      write_xml(xml_data, file=cache_path) 
    }else{
      NULL
    }
  }
  xml_data
}


# get protein locations (e.g. "Cell membrane")
get_locations <- function(uniprot_xml){
  if(is.null(uniprot_xml)) return(NA)
  
  locations <- xml_find_all(uniprot_xml, ".//location")
  locations <- locations[! is.na(xml_attr(locations, "evidence"))]
  unique(xml_text(locations))
}


# get transmembrane regions
get_nr_transmembranes <- function(uniprot_xml){
  if(is.null(uniprot_xml)) return(NA)
  
  features <- xml_find_all(uniprot_xml, ".//feature")
  sel_features <- features[xml_attr(features, "type") == "transmembrane region"]
  if(length(sel_features) >= 3){
    TRUE
  }else{
    FALSE
  }
}


# get glycosylation positions
get_glycosylations <- function(uniprot_xml){
  if(is.null(uniprot_xml)) return(NA)
  
  features <- xml_find_all(uniprot_xml, ".//feature")
  glycosylations <- features[xml_attr(features, "type") == "glycosylation site"]
  glyco_desc <- xml_attr(glycosylations, "description")
  n_linked <- grep("N-linked", glyco_desc)
  
  if(length(n_linked) >= 3){
    TRUE
  }else{
    FALSE
  }
}


# get signal peptide postion
get_signal_pep <- function(uniprot_xml){
  if(is.null(uniprot_xml)) return(NA)
  
  features <- xml_find_all(uniprot_xml, ".//feature")
  signal_pep <- features[xml_attr(features, "type") == "signal peptide"]
  start_end <- xml_attr(xml_children(xml_children(signal_pep)), "position")
  paste(start_end, collapse="-")
}


# get PTM's
get_ptms <- function(uniprot_xml){
  if(is.null(uniprot_xml)) return(NA)
  
  features <- xml_find_all(uniprot_xml, ".//feature")
  ptms <- as.character(features[xml_attr(features, "type") == "modified residue"])
  clean.ptms <- gsub("\\\"|\n|/>", ";", ptms)
  ptm.list <- gsub(".+description=;(.+?);.+", '\\1', clean.ptms, perl=TRUE)
  ptm.unique <- unique(ptm.list)
  
  # remove entries with paranthesis at the start
  ptm.unique <- sub("^\\(.+\\)\\s+", "", ptm.unique)
  
  paste(ptm.unique, collapse=";")
}


# get cross-links (e.g. glycyl lysine)
get_ubiquitin <- function(uniprot_xml, i){
  if(is.null(uniprot_xml)) return(NA)
  
  features <- xml_find_all(uniprot_xml, ".//feature")
  cross_links <- features[xml_attr(features, "type") == "cross-link"]
  cross_link_desc <- unique(xml_attr(cross_links, "description"))
  
  ubiquitins <- sub(".+\\(interchain.+in\\s+(.+)\\).*", "\\1", cross_link_desc)
  all_words <- unlist(strsplit(ubiquitins, "[\\s+|,\\s+]", perl=TRUE)) 
  valid_words <- unique(all_words[all_words %in% valid_ubiquitin])
  
  paste(unique(valid_words), collapse = ";")
}


# get homodimer annotation
get_homodimer <- function(uniprot_xml){
  if(is.null(uniprot_xml)) return(NA)
  
  comments <- xml_find_all(uniprot_xml, ".//comment")
  subunits <- comments[xml_attr(comments, "type") == "subunit"]
  dimers <- grep('homodimer', subunits)
  if(length(dimers) >= 1) TRUE else FALSE
}


# compositionally biased region
get_compositionally_biased_region <- function(uniprot_xml){
  if(is.null(uniprot_xml)) return(NA)
  
  features <- xml_find_all(uniprot_xml, ".//feature")
  cbr <- features[xml_attr(features, "type") == "compositionally biased region"]
  cbr_desc <- xml_attr(cbr, "description")
  cbr_desc_clean <- unique(sub("\\s+\\(.+\\)$", "", perl=TRUE, cbr_desc))
  
  paste(cbr_desc_clean, collapse=";")
}


# see if there is a cleavage with >20AA
get_is_cleaved <- function(uniprot_xml){
  if(is.null(uniprot_xml)) return(NA)
  
  features <- xml_find_all(uniprot_xml, ".//feature")
  chain <- features[xml_attr(features, "type") == "chain"]
  chain_desc <- xml_attr(cbr, "description")
  
  locations <- xml_find_all(chain, ".//location")
  begin_pos <- as.numeric(sub(".+begin position\\=\"(\\d+).*", "\\1", locations))
  end_pos <- as.numeric(sub(".+end position\\=\"(\\d+).*", "\\1", locations))
  
  if(length(chain_desc) > 0) {
    print(chain_desc)
    print(end_pos - begin_pos)
  }
  
  
}

