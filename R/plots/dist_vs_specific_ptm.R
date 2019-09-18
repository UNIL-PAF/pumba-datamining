library(ggplot2)

rm(list=ls())


ptm <- "glycyl"
neg_label <- paste0("no ", ptm)
pos_label <- ptm

# load results
res_path <- ("/Users/admin/tmp/datamining_pumba/results/scop1568714129.32455.RData")
load(res_path)

for(sample in names(results)){
  
  peak_dists <- results[[sample]]
  peak_dists$color <- neg_label
  peak_dists$color[grep(ptm, peak_dists$ptms, ignore.case = TRUE)] <- pos_label
  
  # colors
  col_vec <- c("grey", "red")
  names(col_vec) <- c(neg_label, pos_label)
  
  
  show_mass_range <- c(10, 300)
  
  p <- ggplot(data=peak_dists[order(peak_dists$color),], aes(x=theo_weights, y=closest_peak_masses, colour=color)) 
  p <- p + geom_point(data = peak_dists[peak_dists$color == neg_label,], alpha=0.2)
  p <- p + geom_point(data = peak_dists[peak_dists$color != neg_label,], alpha=0.7)
  p <- p + scale_x_log10()
  p <- p + scale_y_log10()
  p <- p + coord_cartesian(xlim=show_mass_range,ylim=show_mass_range)
  p <- p + geom_abline(intercept=0, slope=1)
  p <- p + scale_colour_manual(values=col_vec)
  p <- p + theme_bw()
  p <- p + ggtitle(sample)
  print(p)
  
  
  # p <- ggplot(data=peak_dists, aes(theo_weights, closest_peak_dists, colour=color)) 
  # p <- p + geom_point(alpha=0.2)
  # # p <- p + geom_text(data=peak_dists[peak_dists$theo_weights > 200, ], aes(theo_weights, closest_peak_dists, label=protein_acs), srt=90, col="red", alpha=0.3)
  # p <- p + xlim(0, 100)
  # p <- p + ylim(-2, 2)
  # p <- p + scale_colour_manual(values=col_vec)
  # p <- p + geom_hline(yintercept = 0)
  # p <- p + theme_bw()
  # p <- p + ggtitle(sample)
  # print(p)
  # 
  # 
  # p <- ggplot(data=peak_dists, aes(color, closest_peak_dists, colour=color)) 
  # p <- p + geom_boxplot()
  # p <- p + geom_hline(yintercept = 0)
  # p <- p + theme_bw()
  # p <- p + ggtitle(sample)
  # print(p)
  
}