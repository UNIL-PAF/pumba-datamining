library(ggplot2)

rm(list=ls())

# load results
res_path <- ("/Users/admin/tmp/datamining_pumba/results/hydrophopicity_1569510778.62979.RData")
load(res_path)
res <- results[['HCT']]

scop_classes <- read.csv("/Users/admin/Work/PAF/projects/pumba/data/datamining/scop/SCOP_classes.csv")

for(scop_class_short in scop_classes$class){
  scop_class_name <- as.character(scop_classes$name[scop_classes$class == scop_class_short])
  
  # Membrane class from SCOP
  res$is_scop_class <- unlist(lapply(res$scop_classes, function(x){ if(length(grep(scop_class_short, x)) > 0) scop_class_name else ""}))
  
  show_mass_range <- c(10, 300)
  
  p <- ggplot(res[order(res$is_scop_class),], aes(x=theo_weights, y=closest_peak_masses, color=is_scop_class)) +
    geom_point(data=res[res$is_scop_class == "", ], alpha=0.1, color="grey") +
    geom_point(data=res[res$is_scop_class != "", ]) +
    scale_x_log10() +
    scale_y_log10() +
    coord_cartesian(xlim=show_mass_range,ylim=show_mass_range) +
    geom_abline(intercept=0, slope=1) +
    ggtitle(scop_class_name) +
    theme_bw()
  
  print(p)
}

