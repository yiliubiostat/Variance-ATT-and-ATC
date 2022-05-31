### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Variance estimations ATE ATT ATC ~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### Sandwich Variance for DR ATT and ATC

### by Yunji Zhou


### Average treatment effect
ATE <- function(y, z, X, DR=FALSE, X.out=NA){
  
  # module for checking balance in the weighted sample
  diff <- function(cov, z, ps){
    # cov: covariate
    # z: treatment status
    # ps: estimated propensity scores
    
    v1 <- cov[z == 1]
    v0 <- cov[z == 0]
    w1 <- 1/ps[z == 1]
    w0 <- 1/(1-ps[z == 0])
    n1 <- length(v1)
    n0 <- length(v0)
    delta <- abs(sum(v1*w1) / sum(w1) - sum(v0*w0) / sum(w0))
    # tstat <- delta / sqrt(var(v1)/n1 + var(v0)/n0)
    tstat <- delta / sqrt(var(v1) + var(v0))
    
    # return the absolute standardized difference
    return(tstat)
  }
  
  # summary statistics
  n1 <- sum(z)             # number of treated
  n0 <- sum(1-z)           # number of untreated
  n <- n0 + n1             # total sample size
  V <- cbind(1, X)         # design matrix (including intercept)
  
  # estimate ps
  fit <- glm(z ~ X, family = binomial(link = "logit"))
  e.h <- as.numeric(fit$fitted.values)
  g.h <- 1
  
  # point estimate - ATE
  mu1.h <- sum(z*y*(1/e.h)) / sum(z*(1/e.h))
  mu0.h <- sum( (1-z)*y* (1/(1-e.h)) ) / sum((1-z)*(1/(1-e.h)))
  tau <- mu1.h - mu0.h
  
  if(DR==FALSE){
  # sandwich variance estimator
  A.11 <- crossprod(sqrt(e.h*(1-e.h)) * V) / n  #(v by v matrix)
  A.22 <- mean(z/e.h*g.h) # scalar 
  A.33 <- mean((1-z)/(1-e.h)*g.h) # scalar
  A.21 <- -apply(z*(0-(1-e.h)*g.h*V)/e.h*(y-mu1.h),2,mean) # vector of length v
  A.31 <- -apply((1-z)*(0+e.h*g.h*V)/(1-e.h)*(y-mu0.h), 2, mean) # vector of length v
  
  A <- rbind(A.11,A.21,A.31)
  A <- cbind(A,c(rep(0,nrow(A.11)),A.22,0),c(rep(0,nrow(A.11)+1),A.33))
  A.inv <- solve(A)
  A.t.inv <- solve(t(A))
  
  psi <- cbind((z-e.h)*V, z*(1/e.h)*(y-mu1.h), (1-z)*(1/(1-e.h))*(y-mu0.h))
  B <- crossprod(psi) / n
  
  sigma <- A.inv %*% B %*% A.t.inv
  c <- c(rep(0,ncol(V)),1,-1)
  
  VAR <- t(c) %*% sigma %*% c /n
  se <- as.numeric(sqrt(VAR))
  
  }else{ 
    
    # Doubly robust estimation
    W <- cbind(1,X.out)
    
    # point estimate
    out.ctrl <- lm(y~.,data = data.frame(y,X.out)[z==0,])
    out.trt  <- lm(y~.,data = data.frame(y,X.out)[z==1,])
    
    m0.h <- predict(out.ctrl,as.data.frame(X.out) )
    m1.h <- predict(out.trt, as.data.frame(X.out) )
    
    w1.h <- (z*g.h/e.h) / sum(z*g.h/e.h)
    w0.h <- ((1-z)*g.h/(1-e.h)) / sum((1-z)*g.h/(1-e.h))
    
    tau <- sum(w1.h*(y-m1.h)) - sum(w0.h*(y-m0.h)) + mean(m1.h-m0.h)  
    
    # sandwich variance estimator
    A.11 <- crossprod(sqrt(e.h*(1-e.h)) * V) / n  #(v by v matrix)
    A.22 <- crossprod(sqrt(z) * W) / n # W by W matrix
    A.33 <- crossprod(sqrt(1-z) * W) / n # W by W matrix
    A.42 <- A.53 <- - apply(W,2,mean) # 1 by w vector
    A.44 <- A.55 <- 1  # scalar or vector?
    A.61 <- ( z*(1-e.h)/e.h*(y-m1.h-sum(w1.h*(y-m1.h))) )%*% V /n # 1 by v vector
    A.62 <- (z/e.h) %*% W /n # 1 by w vector
    A.66 <- mean(z/e.h) # scalar
    A.71 <- -( (1-z)*e.h/(1-e.h)*(y-m0.h-sum(w0.h*(y-m0.h))) ) %*% V /n # 1 by v vector
    A.73 <- ( (1-z)/(1-e.h) ) %*% W / n # 1 by w vector
    A.77 <- mean((1-z)/(1-e.h)) # scalar
    
    
    A.col.1 <- rbind(A.11, matrix(0,nrow = nrow(A.22)+nrow(A.33)+1+1, ncol = ncol(A.11)), A.61, A.71)
    A.col.2 <- rbind(matrix(0,nrow = nrow(A.11),ncol = ncol(A.22)),A.22, matrix(0,nrow = nrow(A.33),ncol = ncol(A.22)), A.42, rep(0, ncol(A.22)), A.62, rep(0, ncol(A.22)))
    A.col.3 <- rbind(matrix(0,nrow = nrow(A.22)+nrow(A.11), ncol = ncol(A.33)), A.33, matrix(0,nrow = 1, ncol = ncol(A.33)), A.53, matrix(0,nrow = 1, ncol = ncol(A.33)), A.73)
    A.col.4 <- c(rep(0, nrow(A.22)+nrow(A.11)+nrow(A.33)),A.44,rep(0,3))
    A.col.5 <- c(rep(0, nrow(A.22)+nrow(A.11)+nrow(A.33)+1),A.55,rep(0,2))
    A.col.6 <- c(rep(0, nrow(A.22)+nrow(A.11)+nrow(A.33)+1+1),A.66,0)
    A.col.7 <- c(rep(0, nrow(A.22)+nrow(A.11)+nrow(A.33)+1+1+1), A.77)
    
    A <- cbind(A.col.1,A.col.2,A.col.3,A.col.4,A.col.5,A.col.6,A.col.7)
    A.inv <- solve(A)
    A.t.inv <- solve(t(A))  
    
    psi <- cbind((z-e.h)*V, c(z  * (y-W %*% out.trt$coefficients)) * W, c( (1-z) * (y-W %*% out.ctrl$coefficients) )* W, 
                 m1.h-mean(m1.h), m0.h-mean(m0.h), z/e.h*(y-m1.h-sum(w1.h*(y-m1.h))), (1-z)/(1-e.h)*(y-m0.h-sum(w0.h*(y-m0.h))))
    B <- crossprod(psi) / n
    
    sigma <- A.inv %*% B %*% A.t.inv
    c <- c(rep(0,ncol(V)+ncol(W)+ncol(W)),1,-1,1,-1)
    
    VAR <- t(c) %*% sigma %*% c /n
    se <- as.numeric(sqrt(VAR))
    
  }
  
  # check balance (expected to be zero)
  asd <- rep(NA, ncol(X))
  for(k in 1:ncol(X)){
    asd[k] <- diff(cov = X[,k], z = z, ps = e.h)
  }
  
  # output the quantities of interest
  return(list(tau=tau, se=se, asd = asd))
}

  
# Simpler version with only point estimator
ATE.PE <- function(y, z, X, DR=FALSE, X.out=NA){
      
  # module for checking balance in the weighted sample
  diff <- function(cov, z, ps){
        # cov: covariate
        # z: treatment status
        # ps: estimated propensity scores
        
        v1 <- cov[z == 1]
        v0 <- cov[z == 0]
        w1 <- 1/ps[z == 1]
        w0 <- 1/(1-ps[z == 0])
        n1 <- length(v1)
        n0 <- length(v0)
        delta <- abs(sum(v1*w1) / sum(w1) - sum(v0*w0) / sum(w0))
        # tstat <- delta / sqrt(var(v1)/n1 + var(v0)/n0)
        tstat <- delta / sqrt(var(v1) + var(v0))
        
        # return the absolute standardized difference
        return(tstat)
      }
      
      # summary statistics
      n1 <- sum(z)             # number of treated
      n0 <- sum(1-z)           # number of untreated
      n <- n0 + n1             # total sample size
      V <- cbind(1, X)         # design matrix (including intercept)
      
      # estimate ps
      fit <- glm(z ~ X, family = binomial(link = "logit"))
      e.h <- as.numeric(fit$fitted.values)
      g.h <- 1
      
      # point estimate - ATE
      mu1.h <- sum(z*y*(1/e.h)) / sum(z*(1/e.h))
      mu0.h <- sum( (1-z)*y* (1/(1-e.h)) ) / sum((1-z)*(1/(1-e.h)))
      tau <- mu1.h - mu0.h
      
      if(DR==FALSE){
        # sandwich variance estimator
        A.11 <- crossprod(sqrt(e.h*(1-e.h)) * V) / n  #(v by v matrix)
        A.22 <- mean(z/e.h*g.h) # scalar 
        A.33 <- mean((1-z)/(1-e.h)*g.h) # scalar
        A.21 <- -apply(z*(0-(1-e.h)*g.h*V)/e.h*(y-mu1.h),2,mean) # vector of length v
        A.31 <- -apply((1-z)*(0+e.h*g.h*V)/(1-e.h)*(y-mu0.h), 2, mean) # vector of length v
        
        A <- rbind(A.11,A.21,A.31)
        A <- cbind(A,c(rep(0,nrow(A.11)),A.22,0),c(rep(0,nrow(A.11)+1),A.33))
        A.inv <- solve(A)
        A.t.inv <- solve(t(A))
        
        psi <- cbind((z-e.h)*V, z*(1/e.h)*(y-mu1.h), (1-z)*(1/(1-e.h))*(y-mu0.h))
        B <- crossprod(psi) / n
        
        sigma <- A.inv %*% B %*% A.t.inv
        c <- c(rep(0,ncol(V)),1,-1)
        
        VAR <- t(c) %*% sigma %*% c /n
        se <- as.numeric(sqrt(VAR))
        
      }else{ 
        
        # Doubly robust estimation
        W <- cbind(1,X.out)
        
        # point estimate
        out.ctrl <- lm(y~.,data = data.frame(y,X.out)[z==0,])
        out.trt  <- lm(y~.,data = data.frame(y,X.out)[z==1,])
        
        m0.h <- predict(out.ctrl,as.data.frame(X.out) )
        m1.h <- predict(out.trt, as.data.frame(X.out) )
        
        w1.h <- (z*g.h/e.h) / sum(z*g.h/e.h)
        w0.h <- ((1-z)*g.h/(1-e.h)) / sum((1-z)*g.h/(1-e.h))
        
        tau <- sum(w1.h*(y-m1.h)) - sum(w0.h*(y-m0.h)) + mean(m1.h-m0.h)  
        se <- NA
        
      }
      
      # check balance (expected to be zero)
      asd <- rep(NA, ncol(X))
      for(k in 1:ncol(X)){
        asd[k] <- diff(cov = X[,k], z = z, ps = e.h)
      }
      
      # output the quantities of interest
      return(list(tau=tau, se=se, asd = asd))
  }

### Average treatment effect on the treated
ATT <- function(y, z, X, DR=FALSE, X.out=NA){
  
  # module for checking balance in the weighted sample
  diff <- function(cov, z, ps){
    # cov: covariate
    # z: treatment status
    # ps: estimated propensity scores
    
    v1 <- cov[z == 1]
    v0 <- cov[z == 0]
    w1 <- rep(1,length(ps[z == 1]))
    w0 <- ps[z == 0]/(1-ps[z == 0])
    n1 <- length(v1)
    n0 <- length(v0)
    delta <- abs(sum(v1*w1) / sum(w1) - sum(v0*w0) / sum(w0))
    # tstat <- delta / sqrt(var(v1)/n1 + var(v0)/n0)
    tstat <- delta / sqrt(var(v1) + var(v0))
    
    # return the absolute standardized difference
    return(tstat)
  }
  
  # summary statistics
  n1 <- sum(z)             # number of treated
  n0 <- sum(1-z)           # number of untreated
  n <- n0 + n1             # total sample size
  V <- cbind(1, X)         # design matrix (including intercept)
  
  # estimate ps
  fit <- glm(z ~ X, family = binomial(link = "logit"))
  e.h <- as.numeric(fit$fitted.values)
  g.h <- e.h
  
  # point estimate - ATT
  mu1.h <- sum(z*y*(1)) / sum(z*(1))
  mu0.h <- sum( (1-z)*y* (e.h/(1-e.h)) ) / sum((1-z)*(e.h/(1-e.h)))
  tau <- mu1.h - mu0.h
  
  if(DR==FALSE){
  # sandwich variance estimator
  A.11 <- crossprod(sqrt(e.h*(1-e.h)) * V) / n  #(v by v matrix)
  A.22 <- mean(z/e.h*g.h) # scalar 
  A.33 <- mean((1-z)/(1-e.h)*g.h) # scalar
  A.21 <- -apply(z*(e.h*(1-e.h)*V-(1-e.h)*g.h*V)/e.h*(y-mu1.h),2,mean) # vector of length v
  A.31 <- -apply((1-z)*(e.h*(1-e.h)*V+e.h*g.h*V)/(1-e.h)*(y-mu0.h), 2, mean) # vector of length v
  
  A <- rbind(A.11,A.21,A.31)
  A <- cbind(A,c(rep(0,nrow(A.11)),A.22,0),c(rep(0,nrow(A.11)+1),A.33))
  A.inv <- solve(A)
  A.t.inv <- solve(t(A))
  
  psi <- cbind((z-e.h)*V, z*1*(y-mu1.h), (1-z)*(e.h/(1-e.h))*(y-mu0.h))
  B <- crossprod(psi) / n
  
  sigma <- A.inv %*% B %*% A.t.inv
  c <- c(rep(0,ncol(V)),1,-1)
  
  VAR <- t(c) %*% sigma %*% c /n
  se <- as.numeric(sqrt(VAR))
  }else{ 
  
  # Doubly robust estimation
  W <- cbind(1,X.out)
  
  # point estimate
  out.ctrl <- lm(y~.,data = data.frame(y,X.out)[z==0,])
  out.trt  <- lm(y~.,data = data.frame(y,X.out)[z==1,])
  
  m0.h <- predict(out.ctrl,as.data.frame(X.out) )
  m1.h <- predict(out.trt, as.data.frame(X.out) )
  
  w1.h <- (z*g.h/e.h) / sum(z*g.h/e.h)
  w0.h <- ((1-z)*g.h/(1-e.h)) / sum((1-z)*g.h/(1-e.h))
  
  tau <- sum( (w1.h-w0.h)*(y-m0.h) )
  
  # sandwich variance estimator
  A.11 <- crossprod(sqrt(e.h*(1-e.h)) * V) / n  #(v by v matrix)
  A.22 <- crossprod(sqrt(1-z) * W) / n # W by W matrix
  A.32 <- z %*% W /n # w by 1 vector
  A.33 <- mean(z) # scalar
  A.41 <- -( (1-z)*e.h/(1-e.h)*(y-m0.h-sum(w0.h*(y-m0.h))) )%*% V /n # v by 1 vector
  A.42 <- ((1-z)*e.h/(1-e.h)) %*% W /n # w by 1 vector
  A.44 <- mean((1-z)*e.h/(1-e.h)) # scalar

  A.col.1 <- rbind(A.11, matrix(0,nrow = nrow(A.22)+1, ncol = ncol(A.11)),A.41)
  A.col.2 <- rbind(matrix(0,nrow = nrow(A.11),ncol = ncol(A.22)),A.22, A.32, A.42)
  A.col.3 <- c(rep(0, nrow(A.22)+nrow(A.11)), A.33, 0)
  A.col.4 <- c(rep(0, nrow(A.22)+nrow(A.11)+1),A.44)
  
  A <- cbind(A.col.1,A.col.2,A.col.3,A.col.4)
  A.inv <- solve(A)
  A.t.inv <- solve(t(A))  
  
  psi <- cbind((z-e.h)*V, c( (1-z)  * (y-W %*% out.ctrl$coefficients)) * W, (z)*(y-m0.h-sum(w1.h*(y-m1.h))), (1-z)*e.h/(1-e.h)*(y-m0.h-sum(w0.h*(y-m0.h))) )
  B <- crossprod(psi) / n
  
  sigma <- A.inv %*% B %*% A.t.inv
  c <- c(rep(0,ncol(V)+ncol(W)),1,-1)
  
  VAR <- t(c) %*% sigma %*% c /n
  se <- as.numeric(sqrt(VAR))
  
  }
  
  # check balance (expected to be zero)
  asd <- rep(NA, ncol(X))
  for(k in 1:ncol(X)){
    asd[k] <- diff(cov = X[,k], z = z, ps = e.h)
  }
  
  # output the quantities of interest
  return(list(tau=tau, se=se, asd = asd))
}


# Simpler version with only point estimator
ATT.PE <- function(y, z, X, DR=FALSE, X.out=NA){
  
  # module for checking balance in the weighted sample
  diff <- function(cov, z, ps){
    # cov: covariate
    # z: treatment status
    # ps: estimated propensity scores
    
    v1 <- cov[z == 1]
    v0 <- cov[z == 0]
    w1 <- rep(1,length(ps[z == 1]))
    w0 <- ps[z == 0]/(1-ps[z == 0])
    n1 <- length(v1)
    n0 <- length(v0)
    delta <- abs(sum(v1*w1) / sum(w1) - sum(v0*w0) / sum(w0))
    # tstat <- delta / sqrt(var(v1)/n1 + var(v0)/n0)
    tstat <- delta / sqrt(var(v1) + var(v0))
    
    # return the absolute standardized difference
    return(tstat)
  }
  
  # summary statistics
  n1 <- sum(z)             # number of treated
  n0 <- sum(1-z)           # number of untreated
  n <- n0 + n1             # total sample size
  V <- cbind(1, X)         # design matrix (including intercept)
  
  # estimate ps
  fit <- glm(z ~ X, family = binomial(link = "logit"))
  e.h <- as.numeric(fit$fitted.values)
  g.h <- e.h
  
  # point estimate - ATT
  mu1.h <- sum(z*y*(1)) / sum(z*(1))
  mu0.h <- sum( (1-z)*y* (e.h/(1-e.h)) ) / sum((1-z)*(e.h/(1-e.h)))
  tau <- mu1.h - mu0.h
  
  if(DR==FALSE){
    # sandwich variance estimator
    A.11 <- crossprod(sqrt(e.h*(1-e.h)) * V) / n  #(v by v matrix)
    A.22 <- mean(z/e.h*g.h) # scalar 
    A.33 <- mean((1-z)/(1-e.h)*g.h) # scalar
    A.21 <- -apply(z*(e.h*(1-e.h)*V-(1-e.h)*g.h*V)/e.h*(y-mu1.h),2,mean) # vector of length v
    A.31 <- -apply((1-z)*(e.h*(1-e.h)*V+e.h*g.h*V)/(1-e.h)*(y-mu0.h), 2, mean) # vector of length v
    
    A <- rbind(A.11,A.21,A.31)
    A <- cbind(A,c(rep(0,nrow(A.11)),A.22,0),c(rep(0,nrow(A.11)+1),A.33))
    A.inv <- solve(A)
    A.t.inv <- solve(t(A))
    
    psi <- cbind((z-e.h)*V, z*1*(y-mu1.h), (1-z)*(e.h/(1-e.h))*(y-mu0.h))
    B <- crossprod(psi) / n
    
    sigma <- A.inv %*% B %*% A.t.inv
    c <- c(rep(0,ncol(V)),1,-1)
    
    VAR <- t(c) %*% sigma %*% c /n
    se <- as.numeric(sqrt(VAR))
  }else{ 
    
    # Doubly robust estimation
    W <- cbind(1,X.out)
    
    # point estimate
    out.ctrl <- lm(y~.,data = data.frame(y,X.out)[z==0,])
    out.trt  <- lm(y~.,data = data.frame(y,X.out)[z==1,])
    
    m0.h <- predict(out.ctrl,as.data.frame(X.out) )
    m1.h <- predict(out.trt, as.data.frame(X.out) )
    
    w1.h <- (z*g.h/e.h) / sum(z*g.h/e.h)
    w0.h <- ((1-z)*g.h/(1-e.h)) / sum((1-z)*g.h/(1-e.h))
    
    tau <- sum( (w1.h-w0.h)*(y-m0.h) )
    
    se <- NA
    
  }
  
  # check balance (expected to be zero)
  asd <- rep(NA, ncol(X))
  for(k in 1:ncol(X)){
    asd[k] <- diff(cov = X[,k], z = z, ps = e.h)
  }
  
  # output the quantities of interest
  return(list(tau=tau, se=se, asd = asd))
}

### Average treatment effect on the control
ATC <- function(y, z, X, DR=FALSE, X.out=NA){
  
  # module for checking balance in the weighted sample
  diff <- function(cov, z, ps){
    # cov: covariate
    # z: treatment status
    # ps: estimated propensity scores
    
    v1 <- cov[z == 1]
    v0 <- cov[z == 0]
    w1 <- (1-ps[z == 1])/ps[z == 1]
    w0 <- rep(1,length(ps[z == 0]))
    n1 <- length(v1)
    n0 <- length(v0)
    delta <- abs(sum(v1*w1) / sum(w1) - sum(v0*w0) / sum(w0))
    # tstat <- delta / sqrt(var(v1)/n1 + var(v0)/n0)
    tstat <- delta / sqrt(var(v1) + var(v0))
    
    # return the absolute standardized difference
    return(tstat)
  }
  
  
  # summary statistics
  n1 <- sum(z)             # number of treated
  n0 <- sum(1-z)           # number of untreated
  n <- n0 + n1             # total sample size
  V <- cbind(1, X)         # design matrix (including intercept)
  
  # estimate ps
  fit <- glm(z ~ X, family = binomial(link = "logit"))
  e.h <- as.numeric(fit$fitted.values)
  g.h <- 1-e.h
  
  # Non-parametric Estimation
  # point estimate
  mu1.h <- sum(z*y*(1-e.h)/e.h) / sum(z*((1-e.h)/e.h))
  mu0.h <- sum( (1-z)*y* (1) ) / sum((1-z)*(1))
  
  
  if(DR==FALSE){
    # Non-parametric Estimation
    tau <- mu1.h - mu0.h
    # sandwich variance estimator
    A.11 <- crossprod(sqrt(e.h*(1-e.h)) * V) / n  #(v by v matrix)
    A.22 <- mean(z/e.h*g.h) # scalar 
    A.33 <- mean((1-z)/(1-e.h)*g.h) # scalar
    A.21 <- -apply(z*(-e.h*(1-e.h)*V-(1-e.h)*g.h*V)/e.h*(y-mu1.h),2,mean) # vector of length v
    A.31 <- -apply((1-z)*(-e.h*(1-e.h)*V+e.h*g.h*V)/(1-e.h)*(y-mu0.h), 2, mean) # vector of length v
    
    A <- rbind(A.11,A.21,A.31)
    A <- cbind(A,c(rep(0,nrow(A.11)),A.22,0),c(rep(0,nrow(A.11)+1),A.33))
    A.inv <- solve(A)
    A.t.inv <- solve(t(A))
    
    psi <- cbind((z-e.h)*V, z*((1-e.h)/e.h)*(y-mu1.h), (1-z)*(1)*(y-mu0.h))
    B <- crossprod(psi) / n
    
    sigma <- A.inv %*% B %*% A.t.inv
    c <- c(rep(0,ncol(V)),1,-1)
    
    VAR <- t(c) %*% sigma %*% c /n
    se <- as.numeric(sqrt(VAR))
    

  }else{ 
    
    # Doubly robust estimation
    W <- cbind(1,X.out)
    
    # point estimate
    out.ctrl <- lm(y~.,data = data.frame(y,X.out)[z==0,])
    out.trt  <- lm(y~.,data = data.frame(y,X.out)[z==1,])
    
    m0.h <- predict(out.ctrl,as.data.frame(X.out) )
    m1.h <- predict(out.trt, as.data.frame(X.out) )
    
    w1.h <- (z*g.h/e.h) / sum(z*g.h/e.h)
    w0.h <- ((1-z)*g.h/(1-e.h)) / sum((1-z)*g.h/(1-e.h))
    
    tau <- sum( (w1.h-w0.h)*(y-m1.h) )
    
    # sandwich variance estimator
    A.11 <- crossprod(sqrt(e.h*(1-e.h)) * V) / n  #(v by v matrix)
    A.22 <- crossprod(sqrt(z) * W) / n # W by W matrix
    A.31 <- ( z*(1-e.h)/e.h*(y-m1.h-sum(w1.h*(y-m1.h))) )%*% V /n # v by 1 vector
    A.32 <- (z*(1-e.h)/e.h) %*% W /n # w by 1 vector
    A.33 <- mean(z*(1-e.h)/e.h) # scalar
    A.42 <- apply((1-z)*W,2,mean) # W by 1 vector
    A.44 <- mean(1-z) # scalar
    
    A.col.1 <- rbind(A.11, matrix(0,nrow = nrow(A.22)+1+1, ncol = ncol(A.11)))
    A.col.2 <- rbind(matrix(0,nrow = nrow(A.11),ncol = ncol(A.22)),A.22, A.32, A.42)
    A.col.3 <- c(rep(0, nrow(A.22)+nrow(A.11)), A.33, 0)
    A.col.4 <- c(rep(0, nrow(A.22)+nrow(A.11)+1),A.44)
    
    A <- cbind(A.col.1,A.col.2,A.col.3,A.col.4)
    A.inv <- solve(A)
    A.t.inv <- solve(t(A))  
    
    psi <- cbind((z-e.h)*V, c(z  * (y-W %*% out.trt$coefficients)) * W, z/e.h*(1-e.h)*(y-m1.h-sum(w1.h*(y-m1.h))), (1-z)*(y-m1.h-sum(w0.h*(y-m0.h))))
    B <- crossprod(psi) / n
    
    sigma <- A.inv %*% B %*% A.t.inv
    c <- c(rep(0,ncol(V)+ncol(W)),1,-1)
    
    VAR <- t(c) %*% sigma %*% c /n
    se <- as.numeric(sqrt(VAR))
    
  }
  
  # check balance (expected to be zero)
  asd <- rep(NA, ncol(X))
  for(k in 1:ncol(X)){
    asd[k] <- diff(cov = X[,k], z = z, ps = e.h)
  }
  
  # output the quantities of interest
  return(list(tau=tau, se=se, asd = asd))
}

# A simpler version with only point estimate
ATC.PE <- function(y, z, X, DR=FALSE, X.out=NA){
  
  # module for checking balance in the weighted sample
  diff <- function(cov, z, ps){
    # cov: covariate
    # z: treatment status
    # ps: estimated propensity scores
    
    v1 <- cov[z == 1]
    v0 <- cov[z == 0]
    w1 <- (1-ps[z == 1])/ps[z == 1]
    w0 <- rep(1,length(ps[z == 0]))
    n1 <- length(v1)
    n0 <- length(v0)
    delta <- abs(sum(v1*w1) / sum(w1) - sum(v0*w0) / sum(w0))
    # tstat <- delta / sqrt(var(v1)/n1 + var(v0)/n0)
    tstat <- delta / sqrt(var(v1) + var(v0))
    
    # return the absolute standardized difference
    return(tstat)
  }
  
  
  # summary statistics
  n1 <- sum(z)             # number of treated
  n0 <- sum(1-z)           # number of untreated
  n <- n0 + n1             # total sample size
  V <- cbind(1, X)         # design matrix (including intercept)
  
  # estimate ps
  fit <- glm(z ~ X, family = binomial(link = "logit"))
  e.h <- as.numeric(fit$fitted.values)
  g.h <- 1-e.h
  
  # Non-parametric Estimation
  # point estimate
  mu1.h <- sum(z*y*(1-e.h)/e.h) / sum(z*((1-e.h)/e.h))
  mu0.h <- sum( (1-z)*y* (1) ) / sum((1-z)*(1))
  
  
  if(DR==FALSE){
    # Non-parametric Estimation
    tau <- mu1.h - mu0.h
    # sandwich variance estimator
    A.11 <- crossprod(sqrt(e.h*(1-e.h)) * V) / n  #(v by v matrix)
    A.22 <- mean(z/e.h*g.h) # scalar 
    A.33 <- mean((1-z)/(1-e.h)*g.h) # scalar
    A.21 <- -apply(z*(-e.h*(1-e.h)*V-(1-e.h)*g.h*V)/e.h*(y-mu1.h),2,mean) # vector of length v
    A.31 <- -apply((1-z)*(-e.h*(1-e.h)*V+e.h*g.h*V)/(1-e.h)*(y-mu0.h), 2, mean) # vector of length v
    
    A <- rbind(A.11,A.21,A.31)
    A <- cbind(A,c(rep(0,nrow(A.11)),A.22,0),c(rep(0,nrow(A.11)+1),A.33))
    A.inv <- solve(A)
    A.t.inv <- solve(t(A))
    
    psi <- cbind((z-e.h)*V, z*((1-e.h)/e.h)*(y-mu1.h), (1-z)*(1)*(y-mu0.h))
    B <- crossprod(psi) / n
    
    sigma <- A.inv %*% B %*% A.t.inv
    c <- c(rep(0,ncol(V)),1,-1)
    
    VAR <- t(c) %*% sigma %*% c /n
    se <- as.numeric(sqrt(VAR))
    
    
  }else{ 
    
    # Doubly robust estimation
    W <- cbind(1,X.out)
    
    # point estimate
    out.ctrl <- lm(y~.,data = data.frame(y,X.out)[z==0,])
    out.trt  <- lm(y~.,data = data.frame(y,X.out)[z==1,])
    
    m0.h <- predict(out.ctrl,as.data.frame(X.out) )
    m1.h <- predict(out.trt, as.data.frame(X.out) )
    
    w1.h <- (z*g.h/e.h) / sum(z*g.h/e.h)
    w0.h <- ((1-z)*g.h/(1-e.h)) / sum((1-z)*g.h/(1-e.h))
    
    tau <- sum( (w1.h-w0.h)*(y-m1.h) )
    se <- NA
    
  }
  
  # check balance (expected to be zero)
  asd <- rep(NA, ncol(X))
  for(k in 1:ncol(X)){
    asd[k] <- diff(cov = X[,k], z = z, ps = e.h)
  }
  
  # output the quantities of interest
  return(list(tau=tau, se=se, asd = asd))
}

