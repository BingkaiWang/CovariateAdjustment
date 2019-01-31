set.seed(100)
library(dplyr)
library(boot)
source("Data_Preprocessing_and_Analysis/adjust_estimator.R")
summary_table_dr <- matrix(NA, nrow = 5, ncol = 3)
rownames(summary_table_dr) <- c("MCI", "METS", "TADS(FLX)", "TADS(CBT)", "TADS(COMB)")
colnames(summary_table_dr) <-c("est", "ci.lower", "ci.upper")

# MCI -----------
load("Data_Preprocessing_and_Analysis/ACDS.rdata")
d <- subset(d, d$arm != "Vitamin E")
y <- d$Y18 - d$Y0
a <- d$arm == "Donepezil"
w <- select(d, age, female, cototscr, mmscore, adtotscr, gdstot, Y0)
delta <- is.na(y)
w <- scale(w)

boot_dr <- boot(data.frame(y, a, w), statistic = function(d, i){
  y <- d[i, 1]
  a <- d[i, 2]
  w <- d[i, -(1:2)]
  delta <- is.na(y)
  adjust_estimator(y, a, w, delta, method = "DR-WLS")
}, R = 5000)
summary_table_dr[1,1] <- boot_dr$t0
summary_table_dr[1,2:3] <- boot.ci(boot_dr, type = "bca")$bca[4:5]

# METS ----------
load("Data_Preprocessing_and_Analysis/METS.rdata")
y <- METS$weightchange
a <- METS$treatment == "Metformin"
w <- select(METS, age, gender, CGI, tobacco, drug, alcohol, weight_baseline, bmi)
w <- scale(w)
delta <- is.na(y)
boot_dr <- boot(data.frame(y, a, w), statistic = function(d, i){
  y <- d[i, 1]
  a <- d[i, 2]
  w <- d[i, -(1:2)]
  delta <- is.na(y)
  adjust_estimator(y, a, w, delta, method = "DR-WLS")
}, R = 5000)
summary_table_dr[2,1] <- boot_dr$t0
summary_table_dr[2,2:3] <- boot.ci(boot_dr, type = "bca")$bca[4:5]

# TADS(FLX) -------
load("Data_Preprocessing_and_Analysis/TADS.rdata")
tads_flx <- filter(tad, treatment %in% c("FLX", "PBO"))
y <- tads_flx$change_score
a <- tads_flx$treatment == "FLX"
w <- select(tads_flx, age, gender, CDRS_baseline,
            CGI, CGAS, RADS, suicide_ideation, depression_episode, comorbidity)
w <- scale(w)
delta <- is.na(y)
boot_dr <- boot(data.frame(y, a, w), statistic = function(d, i){
  y <- d[i, 1]
  a <- d[i, 2]
  w <- d[i, -(1:2)]
  delta <- is.na(y)
  if(sum(delta)){
    adjust_estimator(y, a, w, delta, method = "DR-WLS")
  }else{
    adjust_estimator(y, a, w, method = "ANCOVA")
  }
}, R = 5000)
summary_table_dr[3,1] <- boot_dr$t0
summary_table_dr[3,2:3] <- boot.ci(boot_dr, type = "bca")$bca[4:5]

# TADS(CBT) ------
load("Data_Preprocessing_and_Analysis/TADS.rdata")
tads_cbt <- filter(tad, treatment %in% c("CBT", "PBO"))
y <- tads_cbt$change_score
a <- tads_cbt$treatment == "CBT"
w <- select(tads_cbt, age, gender, CDRS_baseline,
            CGI, CGAS, RADS, suicide_ideation, depression_episode, comorbidity)
w <- scale(w)
delta <- is.na(y)
boot_dr <- boot(data.frame(y, a, w), statistic = function(d, i){
  y <- d[i, 1]
  a <- d[i, 2]
  w <- d[i, -(1:2)]
  delta <- is.na(y)
  adjust_estimator(y, a, w, delta, method = "DR-WLS")
}, R = 5000)
summary_table_dr[4,1] <- boot_dr$t0
summary_table_dr[4,2:3] <- boot.ci(boot_dr, type = "bca")$bca[4:5]

# TADS(COMB)------
load("Data_Preprocessing_and_Analysis/TADS.rdata")
tads_comb <- filter(tad, treatment %in% c("COMB", "PBO"))
y <- tads_comb$change_score
a <- tads_comb$treatment == "COMB"
w <- select(tads_comb, age, gender, CDRS_baseline,
            CGI, CGAS, RADS, suicide_ideation, depression_episode, comorbidity)
w <- scale(w)
delta <- is.na(y)
boot_dr <- boot(data.frame(y, a, w), statistic = function(d, i){
  y <- d[i, 1]
  a <- d[i, 2]
  w <- d[i, -(1:2)]
  delta <- is.na(y)
  adjust_estimator(y, a, w, delta, method = "DR-WLS")
}, R = 5000)
summary_table_dr[5,1] <- boot_dr$t0
summary_table_dr[5,2:3] <- boot.ci(boot_dr, type = "bca")$bca[4:5]

########
round(summary_table_dr, 2)