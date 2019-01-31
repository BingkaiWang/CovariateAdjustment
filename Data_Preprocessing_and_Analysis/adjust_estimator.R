adjust_estimator <- function(y, a, w = NULL, delta = NULL, method = "unadjust", sd.cal = F){
  # Input:
  # y is the outcome, and must be a numerical vector. Here we don't deal with
  # categorical datatype;
  # a is the treatment arm, and must be a binary vector;
  # w is the baseline varible matrix or dataframe;
  # delta is the indicator vector of missing outcomes;
  # method is to specify the estimator to use, it can be "unadjust", "IPW", "DR-WLS", 
  #  "DR-WLS-U" (handling missing value) and "ANCOVA2"
  # Output:
  # And estimator for average treatment effect
  
  
  if(!is.data.frame(w)){
    w <- data.frame(w)
  }
  if(method == "unadjust"){
    est <- mean(y[a == 1]) - mean(y[a == 0])
    if(sd.cal ==  T){
      # sd <- summary(lm(y~a))$coefficients[2,2]
      n <- length(y)
      sd <- sqrt(sum(lm(y~a)$residual^2)/(n-1)/(n-2)/var(a))
    }
  }else if(method == "IPW"){
    propensity <- glm(a~., family = binomial(), data = cbind(a,w))
    g <- predict(propensity, type = "response")
    est <- sum(y * a / g)/sum(a / g) - 
      sum(y * (1-a) / (1-g) )/sum((1-a) / (1-g))
  }else if(method == "DR-WLS"){
    propensity <- glm(arm~., family = binomial(), data = cbind(arm = a,w))
    g <- predict(propensity, type = "response")
    arm1 <- lm(y~., weights = 1/as.vector(g[a == 1]), 
                data = data.frame(y= y[a == 1], w[a == 1, ]))
    arm0 <- lm(y~., weights = 1/(1 - as.vector(g[a == 0])), 
                data = data.frame(y = y[a == 0], w[a == 0, ]))
    est <- mean(predict(arm1, w) - predict(arm0, w))
  }else if(method == "DR-WLS-U"){
    propensity <- glm(arm~., family = binomial(), data = cbind(arm = a,w))
    g <- predict(propensity, type = "response")
    prob_missing <- glm(delta~., family = binomial(), data = cbind(delta, a, w))
    pm <- predict(prob_missing, type = "response")
    indi1 <- (a == 1) & (!is.na(y))
    indi2 <- (a == 0) & (!is.na(y))
    arm1 <- lm(y~., weights = 1/g[indi1]/pm[indi1], 
               data = cbind(y= y[indi1], w[indi1,]))
    arm0 <- lm(y~., weights = 1/(1 - g[indi2])/pm[indi2], 
               data = cbind(y = y[indi2], w[indi2,]))    
    est <- mean(predict(arm1, w) - predict(arm0, w))
  }else if(method == "ANCOVA2"){
    d <- data.frame(y = y, w)
    arm1 <- lm(y~., data = d[a==1,])
    arm0 <- lm(y~., data = d[a==0,])
    est <- mean(predict(arm1, w) - predict(arm0, w))
  }else if(method == "ANCOVA"){
    d <- data.frame(y, a, w)
    est <- lm(y~., data = d)$coef[2]
    if(sd.cal == TRUE){
      # sd <- summary(lm(y~., data = d))$coefficients[2,2]
      n <- length(y)
      p <- ncol(w)
      sd <- sqrt(sum(lm(y~., data = d)$residual^2)/(n-1)/(n-p-2)/var(a))
    }
  }else{
    stop("Undefined method")
  }
  
  if(sd.cal == T){
    return(c(est, sd))
  } else {
    return(est)
  }
}
