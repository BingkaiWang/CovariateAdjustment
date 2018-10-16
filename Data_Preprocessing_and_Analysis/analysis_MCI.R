library(dplyr)
load("Data_Preprocessing_and_Analysis/ACDS.rdata")
source("Data_Preprocessing_and_Analysis/adjust_estimator.R")

d <- subset(d, (d$arm != "Vitamin E") & (!is.na(d$Y18)))
y <- d$Y18 - d$Y0
a <- d$arm == "Donepezil"
w <- select(d, age, female, cototscr, mmscore, adtotscr, gdstot, Y0)
w <- scale(w)
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
MCI_result <- data.frame(treatment_arm = "MCI",
                          aggregate_imbalance = sigma12 %*% solve(sigma22) %*% imbalance/sqrt(sigma12 %*% solve(sigma22) %*% t(sigma12)),
                          unadjusted_estimator = adjust_estimator(y, a),
                          unadj_s.d. = sqrt(sigma11),
                          unadj_c.i. = paste("(", qnorm(0.025,adjust_estimator(y, a),sqrt(sigma11)) %>% round(2), ", ", 
                                             qnorm(0.975,adjust_estimator(y, a),sqrt(sigma11)) %>% round(2), ")", sep = ""),
                          adjusted_estimator = adjust_estimator(y, a, w, method = "ANCOVA"), 
                          adj_s.d. = sqrt(sigma11 - sigma12 %*% solve(sigma22) %*% t(sigma12)),
                          adj_c.i. =  paste("(", qnorm(0.025,adjust_estimator(y, a, w, method = "ANCOVA"), sqrt(sigma11 - sigma12 %*% solve(sigma22) %*% t(sigma12))) %>% round(2), ", ", 
                                            qnorm(0.975,adjust_estimator(y, a, w, method = "ANCOVA"),sqrt(sigma11 - sigma12 %*% solve(sigma22) %*% t(sigma12))) %>% round(2), ")", sep = ""),
                          condi_bias = sigma12 %*% solve(sigma22) %*% imbalance,
                          variance_reduction = sigma12 %*% solve(sigma22) %*% t(sigma12)/sigma11,
                          RE = sigma11/(sigma11 - sigma12 %*% solve(sigma22) %*% t(sigma12))
)
MCI_result
write.csv(MCI_result, file = "DataResults/MCI_result.csv")
