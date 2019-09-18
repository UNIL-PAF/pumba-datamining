library(ggplot2)

rm(list=ls())

# load results
res_path <- ("/Users/admin/tmp/datamining_pumba/results/res_1565182764.29958.RData")
load(res_path)

sample <- "HCT"
res <- results[[sample]]


# pI
res$pI_status <- "neutral"
res$pI_status[res$pIs > 9] <- "pos"
res$pI_status[res$pIs < 5] <- "neg"

show_mass_range <- c(10, 300)

ggplot(res, aes(x=theo_weights, y=closest_peak_masses, color=pI_status)) +
  geom_point(alpha=0.05) +
  geom_point(data=res[res$pI_status == "pos", ]) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim=show_mass_range,ylim=show_mass_range) +
  geom_abline(intercept=0, slope=1) +
  theme_bw()

ggplot(res, aes(x=theo_weights, y=closest_peak_masses, color=pI_status)) +
  geom_point(alpha=0.05) +
  geom_point(data=res[res$pI_status == "neutral", ]) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim=show_mass_range,ylim=show_mass_range) +
  geom_abline(intercept=0, slope=1) +
  theme_bw()



ggplot(res, aes(x=theo_weights, y=closest_peak_masses, color=pI_status)) +
  geom_point(alpha=0.05) +
  geom_point(data=res[res$pI_status == "neg", ]) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim=show_mass_range,ylim=show_mass_range) +
  geom_abline(intercept=0, slope=1) +
  theme_bw()




for(sample in names(results)){
  
  peak_dists <- results[[sample]]
  
  #peak_dists$color <- "neutral"
  #peak_dists$color[peak_dists$charges > 5] <- "pos"
  #peak_dists$color[peak_dists$charges < -5] <- "neg"
  
  peak_dists$color <- "neutral"
  peak_dists$color[peak_dists$pIs > 8] <- "pos"
  peak_dists$color[peak_dists$pIs < 6] <- "neg"
  
  # peak_dists$color <- "neutral"
  # peak_dists$color[peak_dists$charges_by_length > 0.02] <- "pos"
  # peak_dists$color[peak_dists$charges_by_length < -0.02] <- "neg"
  
  p <- ggplot(data=peak_dists, aes(theo_weights, closest_peak_dists, colour=color)) 
  p <- p + geom_point(alpha=0.2)
  #p <- p + geom_text(data=peak_dists[peak_dists$theo_weights > 200, ], aes(theo_weights, closest_peak_dists, label=protein_acs), srt=90, col="red", alpha=0.3)
  p <- p + xlim(0, 100)
  p <- p + ylim(-2, 2)
  p <- p + scale_colour_manual(values=c("pos" = "blue", "neg" = "orange", "neutral" = "grey"))
  p <- p + geom_hline(yintercept = 0)
  p <- p + theme_bw()
  p <- p + ggtitle(sample)
  print(p)
  
  # zoom
  p <- ggplot(data=peak_dists, aes(theo_weights, closest_peak_dists, colour=color)) 
  p <- p + geom_point(alpha=0.2)
  #p <- p + geom_text(data=peak_dists[peak_dists$theo_weights > 200, ], aes(theo_weights, closest_peak_dists, label=protein_acs), srt=90, col="red", alpha=0.3)
  p <- p + xlim(25, 75)
  p <- p + ylim(-2, 2)
  p <- p + scale_colour_manual(values=c("pos" = "blue", "neg" = "orange", "neutral" = "grey"))
  p <- p + geom_hline(yintercept = 0)
  p <- p + theme_bw()
  p <- p + ggtitle(sample)
  print(p)
  
  
  # keep only in between 25 and 75
  flt_peak_dists <- peak_dists[peak_dists$theo_weights >= 25 & peak_dists$theo_weights <= 75, ]
  
  p <- ggplot(data=flt_peak_dists, aes(color, closest_peak_dists, colour=color)) 
  p <- p + geom_boxplot()
  p <- p + geom_hline(yintercept = 0)
  p <- p + theme_bw()
  p <- p + ggtitle(sample)
  print(p)
  
}