

# get the indexes of intensities corresponding to a peak (see "config.R" for thresholds)
get_peak_indexes <- function(ints){
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
      if(! passed_by_zero && ints[i-1] <= peak_detection_zero_thresh) passed_by_zero <- TRUE
    }else{
      if(is_increasing){
        if(ints[i-1] > peak_detection_int_thresh){
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
  
  # filter out peaks below certain threshold relative to strongest peak
  peaks_ints <- ints[peaks_idx]
  min_thresh <- max(peaks_ints) * peak_detection_min_rel_to_strongest
  peaks_ints_ok <- peaks_ints >= min_thresh
  peaks_idx[peaks_ints_ok]
  
}