

# get the indexes of intensities corresponding to a peak (see "config.R" for thresholds)
get_peak_indexes <- function(ints){
  i <- 2
  # test if we're in a increasing part of the curve
  is_increasing <- FALSE
  # assure that we passed close to 0 before looking the next curve
  passed_by_zero <- TRUE
  # if neighboring peaks don't pass at 0 in between, we take the highest one.
  last_int <- 0
  last_peak_int <- 0
  peaks_idx <- c()
  
  while(i <= length(ints)){
    if(ints[i] > ints[i-1]){
      is_increasing <- TRUE
      current_threshold <- last_peak_int * peak_detection_zero_thresh
      if(! passed_by_zero && ints[i-1] <= current_threshold){
        passed_by_zero <- TRUE
      }
    }else{
      if(is_increasing){
        if(ints[i-1] > peak_detection_int_thresh){
          if(passed_by_zero){
            peaks_idx <- c(peaks_idx, i-1)
            last_peak_int <- ints[i-1]
            passed_by_zero <- FALSE
          }else if(ints[i-1] > last_int){
            if(ints[i-1] > ints[peaks_idx[length(peaks_idx)]]){
              last_peak_int <- ints[i-1]
              peaks_idx[length(peaks_idx)] <- i-1 
            }
          }
          last_int <- ints[i-1]
        }
        is_increasing <- FALSE
      }
    }
    i <- i+1
  }
  
  # filter out peaks below certain threshold relative to strongest peak
  peaks_ints <- ints[peaks_idx]

  if(all(is.na(peaks_ints))){
    res <- NULL
  }else{
    min_thresh <- max(peaks_ints, na.rm = TRUE) * peak_detection_min_rel_to_strongest
    peaks_ints_ok <- peaks_ints >= min_thresh
    res <- peaks_idx[peaks_ints_ok]
  }
 
  res
}


get_peak_limits <- function(ints, idx){
  start_int <- ints[idx]
  
  # find the lower end
  current_int <- start_int
  current_idx <- idx
  
  while(current_int > ints[current_idx - 1]){
    current_idx <- current_idx - 1
    current_int <- ints[current_idx]
    if(current_idx == 1) break
  }
  lower_idx <- current_idx
  
  # find the uper end
  current_int <- start_int
  current_idx <- idx
  while(current_int > ints[current_idx + 1]){
    current_idx <- current_idx + 1
    current_int <- ints[current_idx]
    if(current_idx == length(ints)) break
  }
  upper_idx <- current_idx
  
  c(lower_idx, upper_idx)
}