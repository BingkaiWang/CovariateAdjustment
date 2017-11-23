library(Hmisc)
library(dplyr)
aggregate_imbalance <- function(y, a, w){
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
  result <- rbind(imb = imbalance, 
                  proj_coef = sigma12 %*% solve(sigma22),
                  aggr_imb = imbalance * (sigma12 %*% solve(sigma22)) #/as.vector(sqrt(sigma12 %*% solve(sigma22) %*% t(sigma12)))
                  )
  result <- cbind(result, aggregate = c(NA, NA, sum(result[3,])))
  return(result)
}

load("Data_Preprocessing_and_Analysis/TADS.rdata")
d <- subset(tad, !is.na(tad$change_score) & (tad$treatment %in% c("FLX", "PBO")))
TADS_result <- aggregate_imbalance(y = d$change_score, a = d$treatment == "FLX", 
                    w = select(d, age, gender, CDRS_baseline, CGI, CGAS, RADS, suicide_ideation, depression_episode, comorbidity))

load("Data_Preprocessing_and_Analysis/ACDS.rdata")
d <- subset(d, !is.na(d$Y18) & (d$arm %in% c("Placebo", "Donepezil")))
MCI_result <- aggregate_imbalance(y = d$Y18 - d$Y0, a = d$arm == "Donepezil", 
                                  w = select(d, age, female, cototscr, mmscore, adtotscr, gdstot, Y0))

load("Data_Preprocessing_and_Analysis/METS.rdata")
d <- subset(METS, !is.na(METS$weightchange))
METS_result <- aggregate_imbalance(y = d$weightchange, a = d$treatment == "Metformin",
                                   w = select(d, age, gender, CGI, tobacco, drug, alcohol, weight_baseline, bmi))
round(TADS_result,2)
round(MCI_result,2)
round(METS_result,2)
