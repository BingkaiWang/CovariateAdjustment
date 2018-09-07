library(boot)
library(tidyverse)
source("Data_Preprocessing_and_Analysis/adjust_estimator.R")

#####################################################################
# MCI analysis
load("Data_Preprocessing_and_Analysis/ACDS.rdata")
d <- subset(d, (d$arm != "Vitamin E") )
MCI_boot_stat <- function(d, i){
  dd <- d[i,]
  y <- dd$Y18 - dd$Y0
  a <- dd$arm == "Donepezil"
  w <- select(dd, age, female, cototscr, mmscore, adtotscr, gdstot, Y0)
  w <- scale(w)
  delta <- !is.na(y)
  adjust_estimator(y, a, w, delta, method = "DR-WLS-U")
}
boot.out.MCI <- boot(d,MCI_boot_stat,R=1000)
boot.out.MCI$t0
sd(boot.out.MCI$t)
boot.ci(boot.out.MCI,index=1,type="bca")$bca

# METS analysis
load("Data_Preprocessing_and_Analysis/METS.rdata")
METS_boot_stat <- function(d, i){
  dd <- d[i,]
  y <- dd$weightchange
  a <- dd$treatment == "Metformin"
  w <- select(dd, age, gender, CGI, tobacco, drug, alcohol, weight_baseline, bmi)
  w <- scale(w)
  delta <- !is.na(y)
  adjust_estimator(y, a, w, delta, method = "DR-WLS-U")
}
boot.out.METS <- boot(METS,METS_boot_stat,R=1000)
boot.out.METS$t0
boot.ci(boot.out.METS,index=1,type="bca")$bca[4:5]


# TADS
load("Data_Preprocessing_and_Analysis/TADS.rdata")
TADS_boot_stat <- function(data, i, t){
  dd <- data[i, ]
  y <- dd$change_score
  a <- dd$treatment == t
  w <- select(dd, age, gender, CDRS_baseline,
              CGI, CGAS, RADS, suicide_ideation, depression_episode, comorbidity)
  w <- scale(w)
  delta <- !is.na(y)
  adjust_estimator(y, a, w, delta, method = "DR-WLS-U")
}
boot.out.t1 <- boot(tad[tad$treatment %in% c("FLX", "PBO"),],TADS_boot_stat, t = "FLX", R=1000)
boot.out.t1$t0
boot.ci(boot.out.t1,index=1,type="bca")$bca[4:5]

boot.out.t2 <- boot(tad[tad$treatment %in% c("CBT", "PBO"),],TADS_boot_stat, t = "CBT", R=1000)
boot.out.t2$t0
boot.ci(boot.out.t2,index=1,type="bca")$bca[4:5]

boot.out.t3 <- boot(tad[tad$treatment %in% c("COMB", "PBO"),],TADS_boot_stat, t = "COMB", R=1000)
boot.out.t3$t0
boot.ci(boot.out.t3,index=1,type="bca")$bca[4:5]

