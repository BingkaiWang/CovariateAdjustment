library(ggplot2)
library(dplyr)
library(cowplot)
source(file.path("Data_Preprocessing_and_Analysis", "adjust_estimator.R"))

#' Title
#'
#' @param n number of subjects in each simulation
#' @param sim_size number simulations
#' @param num_bins number of bins in the final plot
#' @param range max and min value of aggregate imbalance considered in the final plot
#' @param generating_function data generating function
#' @param plot "scatter" plot or crossbar plot if not specified.
#'
#' @return
#' @export
#'
#' @examples
simulation <- function(n = 200, sim_size = 10000, num_bins = 16, range,
                       generating_function, plot = "crossbar"){
  imbalance <- rep(NA,  sim_size)
  unadj_est <- rep(NA, sim_size)
  adj_est <- rep(NA, sim_size)
  for (m in 1 : sim_size){
    w <- rnorm(n)
    a <- rbinom(n,1,0.5)
    y <- generating_function(a, w) + rnorm(n,0, 0.25)
    imbalance[m] <- mean(w[a == 1]) - mean(w[a == 0])
    unadj_est[m] <- adjust_estimator(y, a, w, method = "unadjust")
    adj_est[m] <- adjust_estimator(y, a, w, method = "ANCOVA")
  }
  if(plot == "scatter"){
    d1 <- data.frame(imbalance = imbalance[1:3000], unadj = unadj_est[1:3000])
    p1 <- ggplot(d1, aes(x = sqrt(n) * imbalance, y = sqrt(n) * unadj)) + 
      geom_point(alpha = 0.5) +
      geom_smooth(method = 'lm',se = F, size = 1.2) + 
      labs(x = expression(paste(sqrt(n), "(imbalance)")),
           y = expression(paste(sqrt(n), "(unadjusted estimator)"))) +
      ylim(-2, 2) + xlim(-5, 5) + 
      theme(text = element_text(size = 20), axis.ticks = element_blank())
    d2 <- data.frame(imbalance = imbalance[1:3000], adj = adj_est[1:3000])
    p2 <- ggplot(d2, aes(x = sqrt(n) * imbalance, y = sqrt(n) * adj)) + 
      geom_point(alpha = 0.5) +
      geom_smooth(method = 'lm',se = F, size = 1.2) + 
      labs(x = expression(paste(sqrt(n), "(imbalance)")),
           y = expression(paste(sqrt(n), "(ANCOVA estimator)"))) +
      ylim(-2, 2) + xlim(-5, 5) + 
      theme(text = element_text(size = 20), axis.ticks = element_blank())    
    scatter_plot <- plot_grid(p1, p2, labels = c("A","B"))
    return(scatter_plot)
  }
  imbalance_bins <- seq(range[1], range[2], by = (range[2]-range[1])/num_bins)
  unadj_cmean = unadj_cvar = adj_cmean = adj_cvar = imbalance_bins
  for (i in 1:length(imbalance_bins)) {
    indi <- imbalance >= imbalance_bins[i] - (range[2]-range[1])/num_bins/2 &
      imbalance < imbalance_bins[i] + (range[2]-range[1])/num_bins/2
    unadj_cmean[i] <- mean(unadj_est[indi])
    unadj_cvar[i] <- var(unadj_est[indi])
    adj_cmean[i] <- mean(adj_est[indi])
    adj_cvar[i] <- var(adj_est[indi])
  }
  summary_unadjust <- cbind(imbalance_bins, unadj_cmean, unadj_cvar)
  summary_adj <- cbind(imbalance_bins, adj_cmean, adj_cvar)
  hist_imbalance <- rbind(summary_unadjust, summary_adj) %>% 
    as.data.frame %>%
    mutate(label = c(rep("unadjusted", nrow(summary_adj)), rep("ANCOVA",nrow(summary_adj))))
  colnames(hist_imbalance) <- c("imbalance", "mean", "variance", "label")
  sim_plot <- ggplot(hist_imbalance, aes(sqrt(n)* imbalance, sqrt(n) *mean, 
                             ymin = sqrt(n) *(mean - sqrt(variance)), 
                             ymax = sqrt(n) * (mean + sqrt(variance)))) +
    geom_crossbar(aes(linetype = label), 
                  fatten = 1, size = 1, width = 0.8*sqrt(n)*(range[2]-range[1])/num_bins, alpha = 1) +
    labs(colour = NULL, 
         x = expression(paste(sqrt(n), "(imbalance)")),
         y = expression(paste(sqrt(n), "(estimator)"))) +
    theme(text = element_text(size = 20), 
          legend.text = element_text(size = 20),
          axis.ticks = element_blank())
  return(sim_plot)
}

p1 <- simulation(sim_size = 10000, num_bins = 16, range = c(-0.3,0.3),
                 generating_function = function(a,w) {0.2*w}) + ylim(-2, 2) + xlim(-5, 5)
pscatter <- simulation(sim_size = 10000, num_bins = 16, range = c(-0.3,0.3),
           generating_function = function(a,w) {0.2*w}, plot = "scatter")
pbarrange <- plot_grid(p1, NULL, labels = c('C', ''), ncol = 2, rel_widths = c(2.62, 1))
plot1 <- plot_grid(pscatter, pbarrange, ncol = 1)
# plot1 <- plot_grid(pscatter, p1, labels = c('', 'C'), ncol = 1, rel_heights = c(1.1,1))
plot1
save_plot('sim1.png', plot1, ncol = 2, nrow = 2)

p2 <- simulation(sim_size = 10000, num_bins = 16, range = c(-0.3,0.3),
           generating_function = function(a,w) {a*w})
p3 <- simulation(sim_size = 10000, num_bins = 16, range = c(-0.3,0.3),
                 generating_function = function(a,w) {0})
p4 <- simulation(sim_size = 10000, num_bins = 16, range = c(-0.3,0.3),
                 generating_function = function(a,w) {(2*a-1)*w})
plot2 <- plot_grid(p1 + theme(plot.title = element_text(size=20),legend.position = c(0,0.9), legend.title = element_blank(), legend.text = element_text(size = 18))+ ggtitle(expression(paste("Y = 0.2W + ", epsilon))), 
                   p3 + theme(plot.title = element_text(size=20),legend.position = "none") + ggtitle(expression(paste("Y = ", epsilon))), 
                   p2 + theme(plot.title = element_text(size=20),legend.position = "none") + ggtitle(expression(paste("Y = AW + ", epsilon))), 
                   p4 + theme(plot.title = element_text(size=20),legend.position = "none") + ggtitle(expression(paste("Y = (2A-1)W + ", epsilon))), 
                   labels = c('1', '2', '3', '4'), ncol = 2, align = 'v', axis = 'l')
plot2
save_plot('sim2.png', plot2, ncol = 2, base_height = 6, base_aspect_ratio = NULL, base_width = 4)
