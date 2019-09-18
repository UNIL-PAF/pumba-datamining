library(ggplot2)

rm(list=ls())

# load results
res_path <- ("/Users/admin/tmp/datamining_pumba/results/frag_1567075972.89237.RData")
load(res_path)

sample <- "HCT"
res <- results[[sample]]

show_mass_range <- c(10, 300)

ggplot(res[order(res$prot_intensities, decreasing=FALSE),], aes(x=theo_weights, y=closest_peak_masses, color=prot_intensities)) +
  geom_point() +
  scale_colour_gradient(low = "white", high = "black", trans="log") +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = show_mass_range, ylim = show_mass_range) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()



# plot intenstiy vs percentage of slices with value
ggplot(res,aes(x=prot_intensities, y=perc_slices)) + 
  geom_point(alpha=0.2, position="jitter")
  #xlim(c(0, 1e10)) 
  #ylim(c(0,0.1))


ggplot(res[order(res$perc_slices, decreasing=FALSE),], aes(x=theo_weights, y=closest_peak_masses, color=perc_slices)) +
  geom_point() +
  scale_colour_gradient(low = "white", high = "black", trans="log") +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = show_mass_range, ylim = show_mass_range) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()


# plot intenstiy vs nr of peaks
ggplot(res,aes(x=prot_intensities, y=nr_peaks)) + 
  geom_point(alpha=0.2, position="jitter")


ggplot(res[order(res$nr_peaks, decreasing=FALSE),], aes(x=theo_weights, y=closest_peak_masses, color=nr_peaks)) +
  geom_point() +
  scale_colour_gradient(low = "#bdd7e7", high = "#08519c") +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = show_mass_range, ylim = show_mass_range) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()

  

# first we label the entries which don't have a good peak at the theoretical mass
res$good_hit <- TRUE
res$good_hit[abs(res$closest_peak_dists) > 0.1] <- FALSE

ggplot(res[order(res$good_hit),], aes(x=theo_weights, y=closest_peak_masses, color=good_hit)) +
  geom_point() +
  #scale_color_manual(values = c("none" = "lightgrey", "both" = "#984ea3", "lower" = "#4daf4a", "higher" = "#ff7f00")) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = show_mass_range, ylim = show_mass_range) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()




# functions to label the results according to the position of the fragments:
# "lower" means that all fragments are below the theoretical mass
# "higher" means that all fragments are above the theoretical mass
# "both" means that some fragments are above and some below theoretical mass
# "none" means that there is only 1 fragment
label_fragments <- function(perc_dists){
  l <- lapply(strsplit(as.character(perc_dists), ","), as.numeric)

  as.factor(unlist(lapply(l, function(dists){
    closest_idx <- which(abs(dists) == min(abs(dists)))
    nr_dists <- length(dists)
    if(nr_dists == 1){
      "none"
    }else if(closest_idx == nr_dists){
      "lower"
    }else if(closest_idx == 1){
      "higher"
    }else{
      "both"
    }
  })))
}

res$fragments <- label_fragments(res$perc_dists)


fragment_order <- c("none", "lower", "higher", "both")

ggplot(res[order(match(res$fragments, fragment_order)),], aes(x=theo_weights, y=closest_peak_masses, color=fragments)) +
  geom_point() +
  scale_color_manual(values = c("none" = "lightgrey", "both" = "#984ea3", "lower" = "#4daf4a", "higher" = "#ff7f00")) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = show_mass_range, ylim = show_mass_range) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()


# get the good hits with lower fragments
lower_frags <- res[res$fragments == "lower" & res$good_hit == TRUE, ]
lower_frags <- lower_frags[order(lower_frags$nr_peaks), ]
lower_frags[,c(3, 1, 12,13,14)]

# get the good hits with higher fragments
higher_frags <- res[res$fragments == "higher" & res$good_hit == TRUE, ]
higher_frags <- higher_frags[order(higher_frags$nr_peaks), ]
higher_frags[,c(3, 1, 12,13,14)]

# check if higher fragment hits could correspond to multimers
find_multimers <- function(peak_masses){
  #multimer_names <- c("dimer", "trimer", "tetramer", "pentamer", "hexamer", "heptamer", "octamer", "nonamer", "decamer", "undecamer", "dodecamer", "tridecamer", "tetradecamer", "pentadecamer", "hexadecamer", "heptadecamer", "octadecamer", "nonadecamer")
  
  l <- lapply(strsplit(as.character(peak_masses), ","), as.numeric)
  
  # set a tolerance of 10%
  multimer_tol <- 0.1
  
  unlist(lapply(l, function(m){
    multimer_res <- c()
    
    if(length(m) > 1){
      for(i in 1:(length(m) - 1)){
        for(k in (i+1):length(m)){
            mass_ratio <- m[k] / m[i]
            rounded_mass_ratio <- round(mass_ratio)
            if(abs(1 - rounded_mass_ratio/mass_ratio) <= multimer_tol & rounded_mass_ratio > 1){
              multimer_res <- c(multimer_res, paste0(rounded_mass_ratio, "-mer"))
            }
        }
      } 
    }
    
    paste(multimer_res, collapse = ",")
  }))
}

res$multimer <- find_multimers(res$peak_masses)






