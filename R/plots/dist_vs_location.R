library(ggplot2)

rm(list=ls())

# location term
location_term <- "membrane"

# load results
res_path <- ("/Users/admin/tmp/datamining_pumba/results/test_1565181053.44648.RData")
load(res_path)

for(sample in names(results)){
  
  peak_dists <- results[[sample]]
  
  neg_location_term <- paste0("not ", location_term)
  peak_dists$color <- neg_location_term
  peak_dists$color[grep(location_term, peak_dists$location, ignore.case = TRUE)] <- location_term

  # colors
  col_vec <- c("grey", "red")
  names(col_vec) <- c(neg_location_term, location_term)
  
  p <- ggplot(data=peak_dists, aes(theo_weights, closest_peak_dists, colour=color)) 
  p <- p + geom_point(alpha=0.2)
  # p <- p + geom_text(data=peak_dists[peak_dists$theo_weights > 200, ], aes(theo_weights, closest_peak_dists, label=protein_acs), srt=90, col="red", alpha=0.3)
  p <- p + xlim(0, 100)
  p <- p + ylim(-2, 2)
  p <- p + scale_colour_manual(values=col_vec)
  p <- p + geom_hline(yintercept = 0)
  p <- p + theme_bw()
  p <- p + ggtitle(sample)
  print(p)
  
}