library(ggplot2)

rm(list=ls())

neg_label <- "low int"
pos_label <- "high int"

# load results
res_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/new_peakdetection_1569505896.27378.RData")
load(res_path)


sample <- "HCT"
res <- results[[sample]]

show_mass_range <- c(10, 300)

ggplot(res[order(res$prot_intensities, decreasing=FALSE),], aes(x=theo_weights, y=closest_peak_masses, color=prot_intensities)) +
  geom_point() +
  scale_colour_gradient(low = "white", high = "#08519c", trans="log") +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = show_mass_range, ylim = show_mass_range) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()

ggplot(res[order(res$prot_intensities, decreasing=FALSE),], aes(x=theo_weights, y=highest_peak_masses, color=prot_intensities)) +
  geom_point() +
  scale_colour_gradient(low = "white", high = "#08519c", trans="log") +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = show_mass_range, ylim = show_mass_range) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()


for(sample in names(results)){
  
  peak_dists <- results[[sample]]
  peak_dists$color <- neg_label
  
  # highlight points >= 3rd quantile
  peak_dists$color[peak_dists$prot_intensities >= summary(peak_dists$prot_intensities)[5]] <- pos_label
  
  # colors
  col_vec <- c("grey", "red")
  names(col_vec) <- c(neg_label, pos_label)
  
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
  
}
  