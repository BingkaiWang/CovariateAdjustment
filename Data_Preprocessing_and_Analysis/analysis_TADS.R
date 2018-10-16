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
  print(diag(cov(w)))
  y1 <- y[a == 1]
  y0 <- y[a == 0]
  w1 <- w[a == 1,]
  w0 <- w[a == 0,]
  n1 <- length(y1)
  n0 <- length(y0)
  imbalance <- colMeans(w1) - colMeans(w0)
  sigma11 <- var(y1)/n1 + var(y0)/n0
  sigma12 <- cov(y1, w1)/n1 + cov(y0, w0)/n0
  sigma22 <- cov(w1)/n1 + cov(w0)/n0
  result <- data.frame(treatment_arm = treatment_arm[1],
                       aggregate_imbalance = sigma12 %*% solve(sigma22) %*% imbalance,
                       unadjusted_estimator = adjust_estimator(y, a),
                       unadj_s.d. = sqrt(sigma11),
                       unadj_c.i. = paste("(", qnorm(0.025,adjust_estimator(y, a),sqrt(sigma11)) %>% round(2), ", ", 
                                          qnorm(0.975,adjust_estimator(y, a),sqrt(sigma11)) %>% round(2), ")", sep = ""),
                       adjusted_estimator = adjust_estimator(y, a, w, method = method), 
                       adj_s.d. = sqrt(sigma11 - sigma12 %*% solve(sigma22) %*% t(sigma12)),
                       adj_c.i. =  paste("(", qnorm(0.025,adjust_estimator(y, a, w, method = method), sqrt(sigma11 - sigma12 %*% solve(sigma22) %*% t(sigma12))) %>% round(2), ", ", 
                                         qnorm(0.975,adjust_estimator(y, a, w, method = method),sqrt(sigma11 - sigma12 %*% solve(sigma22) %*% t(sigma12))) %>% round(2), ")", sep = ""),
                       condi_bias = sigma12 %*% solve(sigma22) %*% imbalance,
                       variance_reduction = sigma12 %*% solve(sigma22) %*% t(sigma12)/sigma11,
                       RE = sigma11/(sigma11 - sigma12 %*% solve(sigma22) %*% t(sigma12))
  )
  result
}

# Complete case analysis of FLX arm to placebo arm
FLX_result <- anaylsis_tad(data = subset(tad, !is.na(tad$CDRS_12)), 
                           treatment_arm = c("FLX", "PBO"), method = "ANCOVA")

# Complete case analysis of CBT arm to placebo arm
CBT_result <- anaylsis_tad(data = subset(tad, !is.na(tad$CDRS_12)), 
                           treatment_arm = c("CBT", "PBO"), method = "ANCOVA")

# Complete case analysis of COMB arm to placebo arm
COMB_result <- anaylsis_tad(data = subset(tad, !is.na(tad$CDRS_12)), 
                           treatment_arm = c("COMB", "PBO"), method = "ANCOVA")

TADS_result <- rbind(FLX_result, CBT_result, COMB_result)
write.csv(TADS_result, file = "DataResults/TADS_result.csv")
