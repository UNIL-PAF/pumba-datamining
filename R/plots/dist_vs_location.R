library(ggplot2)

rm(list=ls())

# location term
location_term <- "mitochond"

# load results
res_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/hydrophopicity_1569510778.62979.RData")
load(res_path)

for(sample in names(results)){
  
  peak_dists <- results[[sample]]
  
  neg_location_term <- paste0("not ", location_term)
  peak_dists$color <- neg_location_term
  peak_dists$color[grep(location_term, peak_dists$location, ignore.case = TRUE)] <- location_term

  # colors
  col_vec <- c("grey", "red")
  names(col_vec) <- c(neg_location_term, location_term)
  
  
  
  show_mass_range <- c(10, 300)
  
  p <- ggplot(data=peak_dists[order(peak_dists$color),], aes(x=theo_weights, y=closest_peak_masses, colour=color)) 
  p <- p + geom_point(data = peak_dists[peak_dists$color == neg_location_term,], alpha=0.2)
  p <- p + geom_point(data = peak_dists[peak_dists$color == location_term,], alpha=0.7)
  p <- p + scale_x_log10()
  p <- p + scale_y_log10()
  p <- p + coord_cartesian(xlim=show_mass_range,ylim=show_mass_range)
  p <- p + geom_abline(intercept=0, slope=1)
  p <- p + scale_colour_manual(values=col_vec)
  p <- p + theme_bw()
  p <- p + ggtitle(sample)
  print(p)
  
  
  # p <- ggplot(data=peak_dists, aes(color, closest_peak_dists, colour=color)) 
  # p <- p + geom_boxplot()
  # p <- p + geom_hline(yintercept = 0)
  # p <- p + theme_bw()
  # p <- p + ggtitle(sample)
  # print(p)
  
}