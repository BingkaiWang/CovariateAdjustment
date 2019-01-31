library(dplyr)
load("Data_Preprocessing_and_Analysis/TADS.rdata")
source("Data_Preprocessing_and_Analysis/adjust_estimator.R")

anaylsis_tad <- function(data, treatment_arm, method){
  # Input:
  # data - a subset of dataframe tad
  # treatment_arm - length-2 vector indicating the two groups to compare, 
  #            with first element as the treatment arm, i.e. c("FLX", "PBO")
  # method - estimator to use. Possible values: "unadjust", "IPW", "DR-WLS", 
  #           "DR-WLS-U" (handling missing value) and "ANCOVA2"
  # Output:
  # A data frame containing information of 
  # treatment_arm, aggregate_imbalance, unadjusted_estimator, unadj_s.d.,
  # adjusted_estimator, adj_s.d., condi_bias, variance_reduction, and RE
  
  data <- data[data$treatment %in% treatment_arm, ]
  y <- data$change_score
  a <- data$treatment == treatment_arm[1]
  w <- select(data, age, gender, CDRS_baseline,
              CGI, CGAS, RADS, suicide_ideation, depression_episode, comorbidity)
  w <- scale(w)
  p <- ncol(w)
  n <- length(y)
  
  # unadjusted ------
  unadj_est <- adjust_estimator(y, a)
  unadj_s.d. <- adjust_estimator(y, a, sd.cal = T)[2]
  unadj_c.i._upper <- unadj_est + abs(qnorm(0.025)) * unadj_s.d.
  unadj_c.i._lower <- unadj_est - abs(qnorm(0.025)) * unadj_s.d.
  
  # ANCOVA -----
  ancova_est <- adjust_estimator(y, a, w, method = "ANCOVA")
  ancova_s.d. <- adjust_estimator(y, a, w, method = "ANCOVA", sd.cal = T)[2]
  ancova_c.i._upper <- ancova_est + abs(qnorm(0.025)) * ancova_s.d.
  ancova_c.i._lower <- ancova_est - abs(qnorm(0.025)) * ancova_s.d.
  
  #result ------
  result <- data.frame(treatment_arm = treatment_arm[1],
                           unadj_est = round(unadj_est, 2),
                           unadj_s.d. = round(unadj_s.d., 2),
                           unadj_c.i. = paste0("(", round(unadj_c.i._lower, 2), ", ", round(unadj_c.i._upper, 2), ")"),
                           ancova_est = round(ancova_est, 2), 
                           ancova_s.d. = round(ancova_s.d., 2),
                           ancova_c.i. =  paste0("(", round(ancova_c.i._lower, 2), ", ", round(ancova_c.i._upper, 2), ")"),
                           variance_reduction = round(1 - ancova_s.d.^2/unadj_s.d.^2, 2),
                           RE = round(unadj_s.d.^2/ancova_s.d.^2, 2)
  )
  
  conditional_bias <- matrix(NA, nrow = 3, ncol = p + 1)
  colnames(conditional_bias) <- c(colnames(w), "Overall")
  rownames(conditional_bias) <- c("coeff.", "Imbal.", "Correction")
  conditional_bias[1,] <- c(lm(y~., data = data.frame(y,a,w))$coef[-(1:2)], NA)
  conditional_bias[2,] <- c(colMeans(w[a == 1,]) - colMeans(w[a == 0,]), NA)
  conditional_bias[3,] <- conditional_bias[1,] * conditional_bias[2,]
  conditional_bias[3, p+1] <- sum(conditional_bias[3, 1:p])
  list(result = result, conditional_bias = round(conditional_bias, 2))
}

# Complete case analysis of FLX arm to placebo arm
FLX_result <- anaylsis_tad(data = subset(tad, !is.na(tad$CDRS_12)), 
                           treatment_arm = c("FLX", "PBO"), method = "ANCOVA")$result

# Complete case analysis of CBT arm to placebo arm
CBT_result <- anaylsis_tad(data = subset(tad, !is.na(tad$CDRS_12)), 
                           treatment_arm = c("CBT", "PBO"), method = "ANCOVA")$result

# Complete case analysis of COMB arm to placebo arm
COMB_result <- anaylsis_tad(data = subset(tad, !is.na(tad$CDRS_12)), 
                           treatment_arm = c("COMB", "PBO"), method = "ANCOVA")$result

TADS_result <- rbind(FLX_result, CBT_result, COMB_result)
TADS_result
write.csv(TADS_result, file = "DataResults/TADS_result.csv")

anaylsis_tad(data = subset(tad, !is.na(tad$CDRS_12)), 
             treatment_arm = c("FLX", "PBO"), method = "ANCOVA")$conditional_bias %>% xtable

anaylsis_tad(data = subset(tad, !is.na(tad$CDRS_12)), 
             treatment_arm = c("CBT", "PBO"), method = "ANCOVA")$conditional_bias %>% xtable

anaylsis_tad(data = subset(tad, !is.na(tad$CDRS_12)), 
             treatment_arm = c("COMB", "PBO"), method = "ANCOVA")$conditional_bias %>% xtable
