load("Data_Preprocessing_and_Analysis/METS.rdata")
source("Data_Preprocessing_and_Analysis/adjust_estimator.R")
METS <- subset(METS, !is.na(METS$weightchange))
y <- METS$weightchange
a <- METS$treatment == "Metformin"
w <- select(METS, age, gender, CGI, tobacco, drug, alcohol, weight_baseline, bmi)
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

# METS result ------
METS_result <- data.frame(treatment_arm = "METS",
                         unadj_est = round(unadj_est, 2),
                         unadj_s.d. = round(unadj_s.d., 2),
                         unadj_c.i. = paste0("(", round(unadj_c.i._lower, 2), ", ", round(unadj_c.i._upper, 2), ")"),
                         ancova_est = round(ancova_est, 2), 
                         ancova_s.d. = round(ancova_s.d., 2),
                         ancova_c.i. =  paste0("(", round(ancova_c.i._lower, 2), ", ", round(ancova_c.i._upper, 2), ")"),
                         variance_reduction = round(1 - ancova_s.d.^2/unadj_s.d.^2, 2),
                         RE = round(unadj_s.d.^2/ancova_s.d.^2, 2)
)
METS_result
write.csv(METS_result, file = "DataResults/METS_result.csv")


# conditional bias ----------
METS_conditional_bias <- matrix(NA, nrow = 3, ncol = p + 1)
colnames(METS_conditional_bias) <- c(colnames(w), "Overall")
rownames(METS_conditional_bias) <- c("coeff.", "Imbal.", "Correction")
METS_conditional_bias[1,] <- c(lm(y~., data = data.frame(y,a,w))$coef[-(1:2)], NA)
METS_conditional_bias[2,] <- c(colMeans(w[a == 1,]) - colMeans(w[a == 0,]), NA)
METS_conditional_bias[3,] <- METS_conditional_bias[1,] * METS_conditional_bias[2,]
METS_conditional_bias[3, p+1] <- sum(METS_conditional_bias[3, 1:p])
round(METS_conditional_bias, 2)

