library(ggplot2)

rm(list=ls())

# load results
res_path <- ("/Users/admin/tmp/datamining_pumba/results/hydrophopicity_1569510778.62979.RData")
load(res_path)

sample <- "HCT"
res <- results[[sample]]

# pI
res$hydrophobicity <- "neutral"
res$hydrophobicity[res$hydrophobicities > 0.3] <- "pos"
res$hydrophobicity[res$hydrophobicities < -0.3] <- "neg"

show_mass_range <- c(10, 300)

ggplot(res, aes(x=theo_weights, y=closest_peak_masses, color=hydrophobicity)) +
  geom_point(alpha=0.05) +
  geom_point(data=res[res$hydrophobicity == "pos", ]) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim=show_mass_range,ylim=show_mass_range) +
  geom_abline(intercept=0, slope=1) +
  theme_bw()

ggplot(res, aes(x=theo_weights, y=closest_peak_masses, color=hydrophobicity)) +
  geom_point(alpha=0.05) +
  geom_point(data=res[res$hydrophobicity == "neutral", ]) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim=show_mass_range,ylim=show_mass_range) +
  geom_abline(intercept=0, slope=1) +
  theme_bw()



ggplot(res, aes(x=theo_weights, y=closest_peak_masses, color=hydrophobicity)) +
  geom_point(alpha=0.05) +
  geom_point(data=res[res$hydrophobicity == "neg", ]) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim=show_mass_range,ylim=show_mass_range) +
  geom_abline(intercept=0, slope=1) +
  theme_bw()


