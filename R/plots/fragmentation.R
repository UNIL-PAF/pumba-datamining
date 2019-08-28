library(ggplot2)

rm(list=ls())

# load results
res_path <- ("/Users/admin/tmp/datamining_pumba/results/cov_1566987107.46999.RData")
load(res_path)

sample <- "HCT"
peak_dists <- results[[sample]]


# plot intenstiy vs percentage of slices with value
ggplot(peak_dists,aes(x=prot_intensities, y=perc_slices)) + 
  geom_point(alpha=0.2, position="jitter") +
  xlim(c(0, 1e10)) 
  #ylim(c(0,0.1))


# plot intenstiy vs nr of peaks
ggplot(peak_dists,aes(x=prot_intensities, y=nr_peaks)) + 
  geom_point(alpha=0.2, position="jitter")
  

summary(peak_dists$prot_intensities)
length(peak_dists$prot_intensities)
sum(peak_dists$prot_intensities == 0)


# look at cases with 4 peaks
peak_dists$protein_acs[peak_dists$nr_peaks == 4]

# O15031, Q8NHP8

peak_dists[peak_dists$nr_peaks == 4, c(3, 1, 12,13,14)]


peak_dists[peak_dists$nr_peaks == 3, c(3, 1, 12,13,14)]
