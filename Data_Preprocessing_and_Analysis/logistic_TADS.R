library(dplyr)
load("Data_Preprocessing_and_Analysis/TADS.rdata")
tad <- subset(tad, !is.na(tad$binary_CGI_improvement) & (tad$treatment %in% c("FLX", "PBO")))
tad$treatment = tad$treatment == "FLX"
tad <- select(tad, binary_CGI_improvement, treatment, age, gender, CDRS_baseline, CGI, 
              CGAS, RADS, suicide_ideation, depression_episode, comorbidity)
# simulation
n <- nrow(tad)
sim_size <- 10000
p1 <- mean(tad$binary_CGI_improvement[tad$treatment])
p0 <- mean(tad$binary_CGI_improvement[!tad$treatment])
q <- (p1 - p0)/(1 - mean(tad$binary_CGI_improvement))
risk_diff <- data.frame(unadjusted = rep(NA, sim_size), standadized = rep(NA, sim_size))
log_odds <- data.frame(unadjusted = rep(NA, sim_size), standardied = rep(NA, sim_size), 
                       logistic = rep(NA, sim_size))
for(i in 1: sim_size){
  # generate data
  stad <- sample_n(tad, size = n, replace = T)
  stad$treatment <- runif(n) < 0.5
  indi <- (stad$binary_CGI_improvement == 0) & (stad$treatment == 1)
  stad$binary_CGI_improvement[indi] <- as.integer(runif(sum(indi)) < q)
  
  # estimation
  glm_result <- glm(binary_CGI_improvement ~ . , data = stad, family = "binomial")
  log_odds$logistic[i] <- glm_result$coefficients[2]
  pred_p1 <- predict(glm_result, cbind(treatment = T, stad[,-(1:2)]), type = "response") %>% mean
  pred_p0 <- predict(glm_result, cbind(treatment = F, stad[,-(1:2)]), type = "response") %>% mean
  risk_diff$standadized[i] <- pred_p1 - pred_p0
  log_odds$standardied[i] <- log(pred_p1/(1-pred_p1)) - log(pred_p0/(1- pred_p0))
  sp1 <- mean(stad$binary_CGI_improvement[stad$treatment])
  sp0 <- mean(stad$binary_CGI_improvement[!stad$treatment])
  risk_diff$unadjusted[i] = sp1 -sp0
  log_odds$unadjusted[i] = log(sp1/(1-sp1)) - log(sp0/(1-sp0))
}

# generate summary table
risk_result <- matrix(NA, nrow = 3, ncol = 2, 
                      dimnames = list(c("mean", "s.e.", "RE"),
                                      c("unadjusted", "standardized")))
logodds_result <- matrix(NA, nrow = 3, ncol = 3, 
                         dimnames = list(c("mean", "s.e.", "RE"),
                                         c("unadjusted", "standardized", "logistic")))
logodds_result[1,] <- colMeans(log_odds)
logodds_result[2,] <- cov(log_odds) %>% diag %>% sqrt
logodds_result[3,] <- logodds_result[1,]^2/logodds_result[2,]^2/(logodds_result[1,3]^2/logodds_result[2,3]^2)
round(logodds_result, 2)
risk_result[1,] <- colMeans(risk_diff)
risk_result[2,] <- cov(risk_diff) %>% diag %>% sqrt
risk_result[3,] <- risk_result[1,]^2/risk_result[2,]^2/(logodds_result[1,3]^2/logodds_result[2,3]^2)
round(risk_result, 2)
