  set.seed(200)
  library(dplyr)
  load("Data_Preprocessing_and_Analysis/ACDS.rdata")
  
  # get binary outcome
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
  d$diagnosis12 <- d$diagnosis12 < 0 # 1 indicates negative and 0 indicate disease
  
  mci <- subset(d, !is.na(d$diagnosis12) & (d$arm %in% c("Placebo", "Donepezil")))
  mci$arm <- mci$arm == "Donepezil"
  # mci <- select(mci, diagnosis12, arm, age, mmscore, apoe4)
  mci <- select(mci, diagnosis12, arm, age, mmscore, apoe4, female, cototscr, adtotscr, gdstot, Y0)
  
  # simulation
  n <- nrow(mci)
  n <- 10*n
  sim_size <- 10000
  p1 <- mean(mci$diagnosis12[mci$arm])
  p0 <- mean(mci$diagnosis12[!mci$arm])
  q <-  (p1 - p0)/(1 - mean(mci$diagnosis12))
  risk_diff <- data.frame(unadjusted = rep(NA, sim_size), standadized = rep(NA, sim_size))
  log_odds <- data.frame(unadjusted = rep(NA, sim_size), standardied = rep(NA, sim_size), 
                         logistic = rep(NA, sim_size))
  simp1 <- rep(NA, sim_size)
  simp2 <- rep(NA, sim_size)
  for(i in 1: sim_size){
    # generate data
    smci <- sample_n(mci, size = n, replace = T)
    smci$arm <- runif(n) < 0.5
    indi <- (smci$diagnosis12 == 0) & (smci$arm == 1)
    smci$diagnosis12[indi] <- as.integer(runif(sum(indi)) < q)
    
    # estimation
    glm_result <- glm(diagnosis12 ~ . , data = smci, family = "binomial")
    log_odds$logistic[i] <- glm_result$coefficients[2]
    pred_p1 <- predict(glm_result, cbind(arm = T, smci[,-(1:2)]), type = "response") %>% mean
    pred_p0 <- predict(glm_result, cbind(arm = F, smci[,-(1:2)]), type = "response") %>% mean
    # pred_p1 <- predict(glm_result, data.frame(treatment = T, CGI = stad[,6]), type = "response") %>% mean
    # pred_p0 <- predict(glm_result, data.frame(treatment = F, CGI = stad[,6]), type = "response") %>% mean
    risk_diff$standadized[i] <- pred_p1 - pred_p0
    log_odds$standardied[i] <- log(pred_p1/(1-pred_p1)) - log(pred_p0/(1- pred_p0))
    sp1 <- mean(smci$diagnosis12[smci$arm])
    sp0 <- mean(smci$diagnosis12[!smci$arm])
    simp1[i] <- sp1
    simp2[i] <- sp0
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
  risk_result[1,] <- colMeans(risk_diff)
  risk_result[2,] <- cov(risk_diff) %>% diag %>% sqrt
  risk_result[3,] <- risk_result[1,]^2/risk_result[2,]^2/(risk_result[1,1]^2/risk_result[2,1]^2)
  risk_result
  log_odds$unadjusted[is.infinite(log_odds$unadjusted)] <- NA
  log_odds <- na.omit(log_odds)
  logodds_result[1,] <- colMeans(log_odds, na.rm =T)
  logodds_result[2,] <- cov(log_odds) %>% diag %>% sqrt
  logodds_result[3,] <- logodds_result[1,]^2/logodds_result[2,]^2/(logodds_result[1,1]^2/logodds_result[2,1]^2)
  logodds_result
  
