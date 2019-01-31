library(dplyr)
load("Data_Preprocessing_and_Analysis/ACDS.rdata")
source("Data_Preprocessing_and_Analysis/adjust_estimator.R")

d <- subset(d, (d$arm != "Vitamin E") & (!is.na(d$Y18)))
y <- d$Y18 - d$Y0
a <- d$arm == "Donepezil"
w <- select(d, age, female, cototscr, mmscore, adtotscr, gdstot, Y0)
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

# MCI result ------
MCI_result <- data.frame(treatment_arm = "MCI",
                         unadj_est = round(unadj_est, 2),
                         unadj_s.d. = round(unadj_s.d., 2),
                         unadj_c.i. = paste0("(", round(unadj_c.i._lower, 2), ", ", round(unadj_c.i._upper, 2), ")"),
                         ancova_est = round(ancova_est, 2), 
                         ancova_s.d. = round(ancova_s.d., 2),
                         ancova_c.i. =  paste0("(", round(ancova_c.i._lower, 2), ", ", round(ancova_c.i._upper, 2), ")"),
                         variance_reduction = round(1 - ancova_s.d.^2/unadj_s.d.^2, 2),
                         RE = round(unadj_s.d.^2/ancova_s.d.^2, 2)
)
MCI_result
write.csv(MCI_result, file = "DataResults/MCI_result.csv")


# conditional bias ----------
MCI_conditional_bias <- matrix(NA, nrow = 3, ncol = p + 1)
colnames(MCI_conditional_bias) <- c(colnames(w), "Overall")
rownames(MCI_conditional_bias) <- c("coeff.", "Imbal.", "Correction")
MCI_conditional_bias[1,] <- c(lm(y~., data = data.frame(y,a,w))$coef[-(1:2)], NA)
MCI_conditional_bias[2,] <- c(colMeans(w[a == 1,]) - colMeans(w[a == 0,]), NA)
MCI_conditional_bias[3,] <- MCI_conditional_bias[1,] * MCI_conditional_bias[2,]
MCI_conditional_bias[3, p+1] <- sum(MCI_conditional_bias[3, 1:p])
round(MCI_conditional_bias, 2)


