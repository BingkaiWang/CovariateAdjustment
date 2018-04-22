library(dplyr)
load("Data_Preprocessing_and_Analysis/METS.rdata")
source("Data_Preprocessing_and_Analysis/adjust_estimator.R")
METS <- subset(METS, !is.na(METS$weightchange))
y <- METS$weightchange
a <- METS$treatment == "Metformin"
w <- select(METS, age, gender, CGI, tobacco, drug, alcohol, weight_baseline, bmi)
w <- scale(w)
y1 <- y[a == 1]
y0 <- y[a == 0]
w1 <- w[a == 1,]
w0 <- w[a == 0,]
n1 <- length(y1)
n0 <- length(y0)
imbalance <- colMeans(w1) - colMeans(w0)
sigma11 <- (var(y1)/n1 + var(y0)/n0)
sigma12 <- (cov(y1, w1)/n1 + cov(y0, w0)/n0)
sigma22 <- (cov(w1)/n1 + cov(w0)/n0)
METS_result <- data.frame(treatment_arm = "METS",
                          unadjusted_estimator = adjust_estimator(y, a),
                          unadj_s.d. = sqrt(sigma11),
                          unadj_c.i. = paste("(", qnorm(0.025,adjust_estimator(y, a),sqrt(sigma11)) %>% round(2), ", ", 
                                             qnorm(0.975,adjust_estimator(y, a),sqrt(sigma11)) %>% round(2), ")", sep = ""),
                          adjusted_estimator = adjust_estimator(y, a, w, method = "ANCOVA"), 
                          adj_s.d. = sqrt(sigma11 - sigma12 %*% solve(sigma22) %*% t(sigma12)),
                          adj_c.i. =  paste("(", qnorm(0.025,adjust_estimator(y, a, w, method = "ANCOVA"), sqrt(sigma11 - sigma12 %*% solve(sigma22) %*% t(sigma12))) %>% round(2), ", ", 
                                            qnorm(0.975,adjust_estimator(y, a, w, method = "ANCOVA"),sqrt(sigma11 - sigma12 %*% solve(sigma22) %*% t(sigma12))) %>% round(2), ")", sep = ""),
                          condi_bias = lm(y~a+w)$coef[-(1:2)] %*% imbalance,
                          variance_reduction = sigma12 %*% solve(sigma22) %*% t(sigma12)/sigma11,
                          RE = sigma11/(sigma11 - sigma12 %*% solve(sigma22) %*% t(sigma12))
)
METS_result

Conditional_bias_decomposition <- rbind(Imbalance = imbalance,
                                        Coefficient = lm(y~a+w)$coef[-(1:2)], 
                                        Contribution_of_conditonal_bias = lm(y~a+w)$coef[-(1:2)] * imbalance)
round(Conditional_bias_decomposition,2)

write.csv(METS_result, file = "DataResults/METS_result.csv")
