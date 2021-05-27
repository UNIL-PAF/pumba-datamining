library(ggplot2)
library(stringr)

rm(list=ls())



# load results
res_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/hydrophopicity_1569510778.62979.RData")
load(res_path)

res <- results[['HCT']]

res$is_glycosylated <- "no glycosylation"
res$is_glycosylated[res$glycosylations != ""] <- "glycosylated"

show_mass_range <- c(10, 300)

ggplot(res[order(res$is_glycosylated, decreasing=TRUE),], aes(x=theo_weights, y=closest_peak_masses, color=is_glycosylated)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("no glycosylation" = "lightgrey", "glycosylated" = "#ff7f00")) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = show_mass_range, ylim = show_mass_range) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()



no_glycosylation_label <- "no glycosylations"
glycosylation_label <- "glycosylated"


  peak_dists$color <- no_glycosylation_label
  peak_dists$color[peak_dists$glycosylations != ""] <- glycosylation_label
  #peak_dists$color[str_count(peak_dists$glycosylations, ",") > 4] <- glycosylation_label
  
  # colors
  col_vec <- c("grey", "red")
  names(col_vec) <- c(no_glycosylation_label, glycosylation_label)
  
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
  
  
  p <- ggplot(data=peak_dists, aes(color, closest_peak_dists, colour=color)) 
  p <- p + geom_boxplot()
  p <- p + geom_hline(yintercept = 0)
  p <- p + theme_bw()
  p <- p + ggtitle(sample)
  print(p)
  