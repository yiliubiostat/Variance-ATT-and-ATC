### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~ PS Methods RHC Data Application ~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### Functions for wild bootstrap variance estimations for DR ATT ATC

### by Yi Liu

### Inverse probability weights
ATE.boot <- function(y, z, X, X.out=NA, R=1000, RV = "Rad"){
  
  # RV: "Rad" or "Exp", the distribution used for generating random vector ksi
  # --- Rad: Rademacher distribution, P(ksi=-1) = P(ksi=1) = 0.5, mean 0
  # --- Exp: exponential(1), mean 1
  
  # summary statistics
  n1 <- sum(z)             # number of treated
  n0 <- sum(1-z)           # number of untreated
  n <- n0 + n1             # total sample size
  V <- cbind(1, X)         # design matrix (including intercept)
  
  # estimate ps
  fit <- glm(z ~ X, family = binomial(link = "logit"))
  e.h <- as.numeric(fit$fitted.values)
  g.h <- 1
  
  # outcome predicts
  out.ctrl <- lm(y~.-y,data = data.frame(y,X.out)[z==0,])
  out.trt  <- lm(y~.-y,data = data.frame(y,X.out)[z==1,])
  
  m0.h <- predict(out.ctrl,as.data.frame(X.out) )
  m1.h <- predict(out.trt, as.data.frame(X.out) )
  
  # weights
  w1.h <- (z*g.h/e.h) / sum(z*g.h/e.h) 
  w0.h <- ((1-z)*g.h/(1-e.h)) / sum((1-z)*g.h/(1-e.h))
  
  # point estimate
  tau <- sum(w1.h*(y-m1.h)) - sum(w0.h*(y-m0.h)) + mean(m1.h-m0.h) 
  
  # Bootstrap: R = 1000 times as default
  seeds <- round(runif(R, min=R*0.01, max=R*100))
  Delta.h <- c()
  
  for(i in 1:R) {
    # random vectors with mean 0 and variance 1
    set.seed(seeds[i])
    
    if(RV=="Exp") ksi <- rexp(n, 1)
    if(RV=="Rad") {rv <- rbinom(n, size=1, prob=0.5); ksi <- ifelse(rv==0, -1, 1)}
    if(!(RV %in% c("Exp", "Rad"))) error("Choose a random vector distribution in Rad or Exp!")
   
    # compute the estimate of Delta
    Delta.h <- c(Delta.h, 
                 sum(ksi*(z/e.h*(y-m1.h)-(1-z)/(1-e.h)*(y-m0.h)+m0.h-m1.h-tau)/sqrt(n)))
  }
  
  z.up <- qnorm(0.75, mean=0, sd=1)
  z.lw <- qnorm(0.25, mean=0, sd=1)
  
  # Bootstrap variance estimate for the limiting distribution
  sigma.boot <- IQR(Delta.h)/(z.up-z.lw)
  
  # Bootstrap variance estimation for tau
  se.boot <- sigma.boot/sqrt(n)
  
  # confidence region
  tstat <- abs(Delta.h)/sigma.boot
  z.crt <- qnorm(0.975, mean=0, sd=1)

  PE.lwr <- tau - z.crt*se.boot
  PE.upr <- tau + z.crt*se.boot
  
  # Bias-corrected estimator and its CI
  tau.pert <- tau + Delta.h/sqrt(n)
  biasC.Est <- 2*tau - mean(tau.pert)
  
  biasC.lwr <- biasC.Est - z.crt*se.boot
  biasC.upr <- biasC.Est + z.crt*se.boot
  
  # output the quantities of interest
  return(list(PE = tau, Std.Boot = se.boot, PE.lwr = PE.lwr, PE.upr = PE.upr,
              biasC.Est = biasC.Est, biasC.lwr = biasC.lwr, biasC.upr = biasC.upr))
}

### Inverse probability weights for the treated
ATT.boot <- function(y, z, X, X.out=NA, R=1000, RV = "Rad"){
  
  # summary statistics
  n1 <- sum(z)             # number of treated
  n0 <- sum(1-z)           # number of untreated
  n <- n0 + n1             # total sample size
  
  # estimate ps
  fit <- glm(z ~ X, family = binomial(link = "logit"))
  e.h <- as.numeric(fit$fitted.values)
  g.h <- e.h
  
  # outcome
  out.ctrl <- lm(y~.-y,data = data.frame(y,X.out)[z==0,])
  out.trt  <- lm(y~.-y,data = data.frame(y,X.out)[z==1,])
  
  m0.h <- predict(out.ctrl,as.data.frame(X.out) )
  m1.h <- predict(out.trt, as.data.frame(X.out) )
  
  # weights
  w1.h <- (z*g.h/e.h) / sum(z*g.h/e.h) 
  w0.h <- ((1-z)*g.h/(1-e.h)) / sum((1-z)*g.h/(1-e.h))
  
  # point estimate
  tau <- sum( (w1.h-w0.h)*(y-m0.h) )
  
  # Bootstrap: R = 1000 times as default
  seeds <- round(runif(R, min=R*0.01, max=R*100))
  Delta.h <- c()
  
  for(i in 1:R) {
    # random vectors with mean 0 and variance 1
    set.seed(seeds[i])
    
    if(RV=="Exp") ksi <- rexp(n, 1)
    if(RV=="Rad") {rv <- rbinom(n, size=1, prob=0.5); ksi <- ifelse(rv==0, -1, 1)}
    if(!(RV %in% c("Exp", "Rad"))) error("Choose a random vector distribution in Rad or Exp!")
    
    # compute the estimate of Delta
    p <- n1/n
    Delta.h <- c(Delta.h, sum(ksi*(1/p*((z-e.h)/(1-e.h)*(y-m0.h)-z*tau)))/sqrt(n) )
  }
  
  z.up <- qnorm(0.75, mean=0, sd=1)
  z.lw <- qnorm(0.25, mean=0, sd=1)
  
  # Bootstrap variance estimate for the limiting distribution
  sigma.boot <- IQR(Delta.h)/(z.up-z.lw)
  
  # Bootstrap variance estimation for tau
  se.boot <- sigma.boot/sqrt(n)
  
  # confidence region
  tstat <- abs(Delta.h)/sigma.boot
  z.crt <- qnorm(0.975, mean=0, sd=1)
  
  PE.lwr <- tau - z.crt*se.boot
  PE.upr <- tau + z.crt*se.boot
  
  # Bias-corrected estimator and its CI
  tau.pert <- tau + Delta.h/sqrt(n)
  biasC.Est <- 2*tau - mean(tau.pert)
  
  biasC.lwr <- biasC.Est - z.crt*se.boot
  biasC.upr <- biasC.Est + z.crt*se.boot
  
  # output the quantities of interest
  return(list(PE = tau, Std.Boot = se.boot, PE.lwr = PE.lwr, PE.upr = PE.upr,
              biasC.Est = biasC.Est, biasC.lwr = biasC.lwr, biasC.upr = biasC.upr))
}


### Inverse probability weights for the control
ATC.boot <- function(y, z, X, X.out=NA, R=1000, RV = "Rad"){
  
  # summary statistics
  n1 <- sum(z)             # number of treated
  n0 <- sum(1-z)           # number of untreated
  n <- n0 + n1             # total sample size
  
  # estimate ps
  fit <- glm(z ~ X, family = binomial(link = "logit"))
  e.h <- as.numeric(fit$fitted.values)
  g.h <- 1-e.h
  
  # outcome
  out.ctrl <- lm(y~.-y,data = data.frame(y,X.out)[z==0,])
  out.trt  <- lm(y~.-y,data = data.frame(y,X.out)[z==1,])
  
  m0.h <- predict(out.ctrl,as.data.frame(X.out) )
  m1.h <- predict(out.trt, as.data.frame(X.out) )
  
  # weights
  w1.h <- (z*g.h/e.h) / sum(z*g.h/e.h) 
  w0.h <- ((1-z)*g.h/(1-e.h)) / sum((1-z)*g.h/(1-e.h))
  
  # point estimate
  tau <- sum( (w1.h-w0.h)*(y-m1.h) )
  
  # Bootstrap: R = 1000 times as default
  seeds <- round(runif(R, min=R*0.01, max=R*100))
  Delta.h <- c()
  
  for(i in 1:R) {
    # random vectors with mean 0 and variance 1
    set.seed(seeds[i])
    
    if(RV=="Exp") ksi <- rexp(n, 1)
    if(RV=="Rad") {rv <- rbinom(n, size=1, prob=0.5); ksi <- ifelse(rv==0, -1, 1)}
    if(!(RV %in% c("Exp", "Rad"))) error("Choose a random vector distribution in Rad or Exp!")
    
    # compute the estimate of Delta
    p <- n1/n
    Delta.h <- c(Delta.h, sum(ksi*(1/(1-p)*(z-e.h)/e.h*(y-m1.h)-(1-z)*tau)/sqrt(n)))
  }
  
  z.up <- qnorm(0.75, mean=0, sd=1)
  z.lw <- qnorm(0.25, mean=0, sd=1)
  
  # Bootstrap variance estimate for the limiting distribution
  sigma.boot <- IQR(Delta.h)/(z.up-z.lw)
  
  # Bootstrap variance estimation for tau
  se.boot <- sigma.boot/sqrt(n)
  
  # confidence region
  tstat <- abs(Delta.h)/sigma.boot
  z.crt <- qnorm(0.975, mean=0, sd=1)
  
  PE.lwr <- tau - z.crt*se.boot
  PE.upr <- tau + z.crt*se.boot
  
  # Bias-corrected estimator and its CI
  tau.pert <- tau + Delta.h/sqrt(n)
  biasC.Est <- 2*tau - mean(tau.pert)
  
  biasC.lwr <- biasC.Est - z.crt*se.boot
  biasC.upr <- biasC.Est + z.crt*se.boot
  
  # output the quantities of interest
  return(list(PE = tau, Std.Boot = se.boot, PE.lwr = PE.lwr, PE.upr = PE.upr,
              biasC.Est = biasC.Est, biasC.lwr = biasC.lwr, biasC.upr = biasC.upr))
}

