library(ggplot2)
library(dplyr)
library(cowplot)
source(file.path("Data_Preprocessing_and_Analysis", "adjust_estimator.R"))

simulation <- function(n = 200, sim_size = 10000, num_bins = 16, range, 
                       generating_function, plot = "crossbar"){
  imbalance <- rep(NA,  sim_size)
  unadj_est <- rep(NA, sim_size)
  adj_est <- rep(NA, sim_size)
  for (m in 1 : sim_size){
    w <- rnorm(n)
    a <- rbinom(n,1,0.5)
    f <- generating_function(a, w)
    y <- f + rnorm(n,0, 0.01)
    imbalance[m] <- mean(w[a == 1]) - mean(w[a == 0])
    unadj_est[m] <- adjust_estimator(y, a, w, method = "unadjust")
    adj_est[m] <- adjust_estimator(y, a, w, method = "ANCOVA2")
  }
  if(plot == "scatter"){
    d1 <- data.frame(imbalance = imbalance[1:3000], unadj = unadj_est[1:3000])
    p1 <- ggplot(d1, aes(x = imbalance, y = unadj)) + 
      geom_point(alpha = 0.5) +
      geom_smooth(method = 'lm',se = F, size = 1.2) + 
      labs(x = "imbalance", y = "unadjusted estimator")
    d2 <- data.frame(imbalance = imbalance[1:3000], adj = adj_est[1:3000])
    p2 <- ggplot(d2, aes(x = imbalance, y = adj)) + 
      geom_point(alpha = 0.5) +
      geom_smooth(method = 'lm',se = F, size = 1.2) + 
      labs(x = "imbalance", y = "adjusted estimator")      
    scatter_plot <- plot_grid(p1, p2, labels = c("A","B"))
    return(scatter_plot)
  }
  aggr_imb <- imbalance * cov(y,w)
  imbalance_bins <- seq(range[1], range[2], by = (range[2]-range[1])/num_bins)
  unadj_cmean = unadj_cvar = adj_cmean = adj_cvar = imbalance_bins
  for (i in 1:length(imbalance_bins)) {
    indi <- aggr_imb >= imbalance_bins[i] - (range[2]-range[1])/num_bins/2 &
      aggr_imb < imbalance_bins[i] + (range[2]-range[1])/num_bins/2
    unadj_cmean[i] <- mean(unadj_est[indi])
    unadj_cvar[i] <- var(unadj_est[indi])
    adj_cmean[i] <- mean(adj_est[indi])
    adj_cvar[i] <- var(adj_est[indi])
  }
  summary_unadjust <- cbind(imbalance_bins, unadj_cmean, unadj_cvar)
  summary_adj <- cbind(imbalance_bins, adj_cmean, adj_cvar)
  hist_imbalance <- rbind(summary_unadjust, summary_adj) %>% 
    as.data.frame %>%
    mutate(label = c(rep("unadjusted", nrow(summary_adj)), rep("adjusted",nrow(summary_adj))))
  colnames(hist_imbalance) <- c("imbalance", "mean", "variance", "label")
  sim_plot <- ggplot(hist_imbalance, aes(sqrt(n)* imbalance, sqrt(n) *mean, 
                             ymin = sqrt(n) *(mean - sqrt(variance)), 
                             ymax = sqrt(n) * (mean + sqrt(variance)))) +
    geom_crossbar(aes(color = label), #position = position_dodge(0),
                  fatten = 1, size = 1, width = 0.8*sqrt(n)*(range[2]-range[1])/num_bins, alpha = 1) +
    labs(colour = NULL, y = "sqrt(n) * estimator") +
    theme(text = element_text(size = 16), legend.text = element_text(size = 16))
  return(sim_plot)
}

simulation(sim_size = 10000, num_bins = 16, range = c(-0.16,0.16),
           generating_function = function(a,w) {a*w})

simulation(sim_size = 10000, num_bins = 16, range = c(-0.16,0.16),
           generating_function = function(a,w) {a*w}, plot = "scatter")
save_plot('scatter.png', scatter_plot, ncol = 2)

simulation(sim_size = 10000, num_bins = 16, range = c(-0.16,0.16),
           generating_function = function(a,w) {a*w + (1-a)*(1-w^2)})
