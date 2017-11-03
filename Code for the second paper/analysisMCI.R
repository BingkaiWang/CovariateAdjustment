set.seed(100)

# Comparing adjusted, unadjusted, and logistic regression estimator
library(dplyr)
library(ggplot2)
library(Hmisc)
library(ltmle)

load("Data_Preprocessing_and_Analysis/ACDS.rdata") # d

# Compute Outcome
d <- mutate(d, diagnosis12 = NA)
diagnosis <- read.csv("Data_Preprocessing_and_Analysis/diagsum.csv")
diagnosis <- select(diagnosis, ptno, viscode, dgalzhei)
for(i in 1: nrow(d)){
  temp <- filter(diagnosis, ptno == d$rid[i]) %>%
    filter(viscode == 'm12')
  if(nrow(temp)){
    d$diagnosis12[i] <- temp$dgalzhei
  }
}
d$diagnosis12[d$diagnosis12 == -1] <- NA
d$diagnosis12 <- d$diagnosis12 > 0

# remove missing values and extract FLX and placebo arm
d <- subset(d, !is.na(d$diagnosis12) & (d$arm %in% c("Placebo", "Donepezil")))


# Baseline variables
W = select(d, age, mmscore, apoe4, female, cototscr, adtotscr, gdstot, Y0)

# Treatment
A = as.integer(d$arm == "Donepezil")

# Outcome
Y = as.integer(!d$diagnosis12)

# Creating dataframe used
data.used = data.frame(Y, A, W)

# Create a funtion that calculates the estimator as well as a 95% CI using the
# percentile bootstrap for a given dataset as well as number of bootstrap replictes
one.sim = function(data.used, n.boot){
  
  # Calculating treatment arm specific means
  res.1.unad = mean(data.used$Y[data.used$A == 1])
  res.0.unad = mean(data.used$Y[data.used$A == 0])
  
  # Fitting G computation and logistic regression estimator 
  log.reg = glm(Y ~., data = data.used, family = "binomial")
  log.reg.1 = log.reg$coefficients[["A"]]
  var.log.reg = summary(log.reg)$coefficients["A", "Std. Error"]^2
  data.a.1 = data.used
  data.a.1$A = 1
  data.a.0 = data.used
  data.a.0$A = 0
  pred.1 =  predict.glm(log.reg, newdata = data.a.1[, -1], type = "response")
  pred.0 =  predict.glm(log.reg, newdata = data.a.0[, -1], type = "response")
  res.gcomp = mean(pred.1) - mean(pred.0)
  
  #boot.0.tmle = rep(NA, n.boot)
  boot.1.unad = rep(NA, n.boot)
  boot.0.unad = rep(NA, n.boot)
  boot.gcomp = rep(NA, n.boot)
  boot.lr = rep(NA, n.boot)
  
  for(i in 1:n.boot){
    bs = sample(1:nrow(data.used), size = nrow(data.used), replace = TRUE)
    boot.1.unad[i] = mean(data.used[bs, ]$Y[data.used[bs, ]$A == 1])
    boot.0.unad[i] = mean(data.used[bs, ]$Y[data.used[bs, ]$A == 0])
    
    # Fitting G computation and logistic regression estimator 
    log.reg.bs = glm(Y ~., data = data.used[bs, ], family = "binomial")
    boot.lr[i] = log.reg.bs$coefficients[["A"]]
    data.a.1 = data.used[bs, ]
    data.a.1$A = 1
    data.a.0 = data.used[bs, ]
    data.a.0$A = 0
    pred.1.bs =  predict.glm(log.reg.bs, newdata = data.a.1[, -1], type = "response")
    pred.0.bs =  predict.glm(log.reg.bs, newdata = data.a.0[, -1], type = "response")
    boot.gcomp[i] = mean(pred.1.bs) - mean(pred.0.bs)
  }
  ret = list(res.gcomp = res.gcomp, var.gcomp = var(boot.gcomp), res.unad = res.1.unad - res.0.unad, var.unad = var(boot.1.unad - boot.0.unad), res.log.reg = log.reg.1, var.log.reg = var.log.reg, var.log.reg.boot = var(boot.lr))
  return(ret)
}

logit = function(x){
  return(log(x/(1-x)))
}

# Create a funtion that calculates the estimator as well as a 95% CI using the
# percentile bootstrap for a given dataset as well as number of bootstrap replictes
one.sim.log.odds = function(data.used, n.boot){
  
  # Calculating treatment arm specific means
  res.1.unad = logit(mean(data.used$Y[data.used$A == 1]))
  res.0.unad = logit(mean(data.used$Y[data.used$A == 0]))
  
  # Fitting G computation and logistic regression estimator 
  log.reg = glm(Y ~., data = data.used, family = "binomial")
  log.reg.1 = log.reg$coefficients[["A"]]
  var.log.reg = summary(log.reg)$coefficients["A", "Std. Error"]^2
  data.a.1 = data.used
  data.a.1$A = 1
  data.a.0 = data.used
  data.a.0$A = 0
  pred.1 =  predict.glm(log.reg, newdata = data.a.1[, -1], type = "response")
  pred.0 =  predict.glm(log.reg, newdata = data.a.0[, -1], type = "response")
  res.gcomp = logit(mean(pred.1)) - logit(mean(pred.0))
  
  #boot.0.tmle = rep(NA, n.boot)
  boot.1.unad = rep(NA, n.boot)
  boot.0.unad = rep(NA, n.boot)
  boot.gcomp = rep(NA, n.boot)
  boot.lr = rep(NA, n.boot)
  
  for(i in 1:n.boot){
    bs = sample(1:nrow(data.used), size = nrow(data.used), replace = TRUE)
    boot.1.unad[i] = logit(mean(data.used[bs, ]$Y[data.used[bs, ]$A == 1]))
    boot.0.unad[i] = logit(mean(data.used[bs, ]$Y[data.used[bs, ]$A == 0]))
    
    # Fitting G computation and logistic regression estimator 
    log.reg.bs = glm(Y ~., data = data.used[bs, ], family = "binomial")
    boot.lr[i] = log.reg.bs$coefficients[["A"]]
    data.a.1 = data.used[bs, ]
    data.a.1$A = 1
    data.a.0 = data.used[bs, ]
    data.a.0$A = 0
    pred.1.bs =  predict.glm(log.reg.bs, newdata = data.a.1[, -1], type = "response")
    pred.0.bs =  predict.glm(log.reg.bs, newdata = data.a.0[, -1], type = "response")
    boot.gcomp[i] = logit(mean(pred.1.bs)) - logit(mean(pred.0.bs))
  }
  
  ret = list(res.gcomp = res.gcomp, var.gcomp = var(boot.gcomp), res.unad = res.1.unad - res.0.unad, var.unad = var(boot.1.unad - boot.0.unad), res.log.reg = log.reg.1, var.log.reg = var.log.reg, var.log.reg.boot = var(boot.lr))
  return(ret)
}

# Running code on non-
result = one.sim.log.odds(data.used, 5000)

# CI
conf.int.gcomp = c(result$res.gcomp - 1.96 * sqrt(result$var.gcomp), result$res.gcomp + 1.96 * sqrt(result$var.gcomp))
conf.int.unad = c(result$res.unad - 1.96 * sqrt(result$var.unad), result$res.unad + 1.96 * sqrt(result$var.unad))
conf.int.log.reg = c(result$res.log.reg - 1.96 * sqrt(result$var.log.reg.boot), result$res.log.reg + 1.96 * sqrt(result$var.log.reg.boot))
#conf.int.gcomp.2 = c(result$res.gcomp - quantile(sqrt(result$var.gcomp), 0.025), result$res.gcomp + quantile(sqrt(result$var.gcomp), 0.975))
#conf.int.unad.2 = c(result$res.unad - quantile(sqrt(result$var.unad), 0.025), result$res.unad + quantile(sqrt(result$var.unad), 0.025))


# P-value
test.stat.gcomp = result$res.gcomp/sqrt(result$var.gcomp)
p.val.gcomp = 2* (1- pnorm(test.stat.gcomp))

test.stat.unad = result$res.unad/sqrt(result$var.unad)
p.val.unad = 2* (1- pnorm(test.stat.unad))

test.stat.logreg = result$res.log.reg/sqrt(result$var.log.reg.boot)
p.val.log.reg = 2* (1- pnorm(test.stat.logreg))

# summary log odds result
CI_num2str <- function(ci){ paste0("(", round(ci[1],2), ", ", round(ci[2], 2), ")")}
output_logodss <- data.frame(estimate_logodds = round(c(result$res.unad, result$res.gcomp, result$res.log.reg),2),
                             se_logodds = sqrt(c(result$var.unad, result$var.gcomp, result$var.log.reg.boot)),
                             CI_logodds = c(CI_num2str(conf.int.unad), CI_num2str(conf.int.gcomp), CI_num2str(conf.int.log.reg)),
                             Pvalue_logodds = c(p.val.unad, p.val.gcomp, p.val.log.reg))


# Running code on non-
result = one.sim(data.used, 5000)

# CI
conf.int.gcomp = c(result$res.gcomp - 1.96 * sqrt(result$var.gcomp), result$res.gcomp + 1.96 * sqrt(result$var.gcomp))
conf.int.unad = c(result$res.unad - 1.96 * sqrt(result$var.unad), result$res.unad + 1.96 * sqrt(result$var.unad))
conf.int.log.reg = c(result$res.log.reg - 1.96 * sqrt(result$var.log.reg.boot), result$res.log.reg + 1.96 * sqrt(result$var.log.reg.boot))
#conf.int.gcomp.2 = c(result$res.gcomp - quantile(sqrt(result$var.gcomp), 0.025), result$res.gcomp + quantile(sqrt(result$var.gcomp), 0.975))
#conf.int.unad.2 = c(result$res.unad - quantile(sqrt(result$var.unad), 0.025), result$res.unad + quantile(sqrt(result$var.unad), 0.025))


# P-value
test.stat.gcomp = result$res.gcomp/sqrt(result$var.gcomp)
p.val.gcomp = 2* (1- pnorm(test.stat.gcomp))

test.stat.unad = result$res.unad/sqrt(result$var.unad)
p.val.unad = 2* (1- pnorm(test.stat.unad))

test.stat.logreg = result$res.log.reg/sqrt(result$var.log.reg.boot)
p.val.log.reg = 2* (1- pnorm(test.stat.logreg))

# summary riskdiff result
CI_num2str <- function(ci){ paste0("(", round(ci[1],2), ", ", round(ci[2], 2), ")")}
output_riskdiff <- data.frame(estimate_riskdiff = round(c(result$res.unad, result$res.gcomp, result$res.log.reg),2),
                              se_riskdiff = sqrt(c(result$var.unad, result$var.gcomp, result$var.log.reg.boot)),
                              CI_riskdiff = c(CI_num2str(conf.int.unad), CI_num2str(conf.int.gcomp), CI_num2str(conf.int.log.reg)),
                              Pvalue_riskdiff = c(p.val.unad, p.val.gcomp, p.val.log.reg))

# combine summary table of riskdiff and logodds
output <- cbind(output_riskdiff, output_logodss)
rownames(output) <- c("Unadjusted", "Standardized", "Logistic")
latexTabular(output)
output
