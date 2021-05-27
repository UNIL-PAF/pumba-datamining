library(ggplot2)

rm(list=ls())

# load results
res_path <- ("/Users/rmylonas/tmp/datamining_pumba/results/int_1565269689.97661.RData")
load(res_path)

sample <- "HCT"
peak_dists <- results[[sample]]

# order data by intensity and add column indicating if high/low/medium intensity
one_third_nr <- round(nrow(peak_dists)/3)
sorted_ints <- sort(peak_dists$prot_intensities)
first_limit <- sorted_ints[one_third_nr]
second_limit <- sorted_ints[one_third_nr * 2]
peak_dists$int_label <- "medium"
peak_dists$int_label[peak_dists$prot_intensities < first_limit] <- "low"
peak_dists$int_label[peak_dists$prot_intensities >= second_limit] <- "high"

ggplot(peak_dists,aes(x=pIs, color=int_label, fill=int_label)) + 
  geom_histogram()

ggplot(peak_dists,aes(x=pIs, color=int_label, fill=int_label)) + 
  geom_histogram(alpha=0.2, position="identity")

ggplot(peak_dists,aes(x=pIs, color=int_label)) + 
  geom_density()





