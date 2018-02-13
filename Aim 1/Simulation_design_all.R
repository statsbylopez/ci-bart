
#### #### #### #### #### #### #### 
#### Simulation design 
#### #### #### #### #### #### #### 



inv_logit = function(x) exp(x) / (1 + exp(x))

# new data generation model
data_generation_all = function(n, ratio, p=10, Zmodel, Ymodel, Aligned) {
  
  #n <- 100000
  #Zmodel <- "Linear"
  #Ymodel <- "Non-Parallel"
  #Aligned=FALSE
  
  
  # total sample size, n=n1+n2+n3
  # n = {700, 3500}
  
  # approximate ratio of treatment groups
  # T1:T2:T3 = 1:2:4, 1:1:1, 4:2:1 
  # ratio = 2, 1, 0.5
  ratio <- 0.5
  if (ratio == 2){
    alpha1 <- -1.8
    alpha2 <- -1.0}
  if (ratio == 1){
    alpha1 <- 0
    alpha2 <- 0}
  if (ratio == 0.5){
    alpha1 <- 1.8
    alpha2 <- 1}

  
  
  
  # number of confounding factors
  # p = {10,50}
  
  # number of relevant confounding factors, p = 9
  # vary across different scenarios (aligned or not as aligned)
  X = matrix(rnorm(p*n), nrow=n, ncol=p)
  x1 = X[,1]; x2 = X[,2]; x3 = X[,3]
  x7 = X[,7]; x8 = X[,8]; x9 = X[,9]
  
  # Z-Model, treatment assignment
  if (Zmodel=="Linear") {
    # linear
    ex1 = exp(alpha1 + 0.5*x1 + 0.7*x2 + 0.5*x3 + 0.8*x7 + 0.2*x8 + 0.8*x9)
    ex2 = exp(alpha2 + 0.9*x1 + 0.3*x2 + 0.9*x3 + 0.2*x7 + 0.6*x8 + 0.2*x9)
    Zp1 = ex1 / (1 + ex1 + ex2)
    Zp2 = ex2 / (1 + ex1 + ex2)
    Zp3 = 1 - Zp1 - Zp2
    
    Zp1 = ifelse(Zp1 < 10^(-10), 0, Zp1)
    Zp2 = ifelse(Zp2 < 10^(-10), 0, Zp2)
    Zp3 = ifelse(Zp3 < 10^(-10), 0, Zp3)
    
  } else if (Zmodel=="Nonlinear") {
    # Nonlinear
    ex1 = exp(alpha1 + 0.5*x1 + 1.0*x2 + 2.0*x3 + 0.2*x1^2 + 0.5*x2^2 + 1.0*x3^2 + 0.6*x7 + 1.2*x8 + 1.5*x9 + 1.0*x7^2 + 0.6*x8^2 + 0.4*x9^3 + 0.5*x1*x7 + 0.7*x2*x8 + 0.9*x3*x9 + 0.5*x7*x8*x9)
    ex2 = exp(alpha2 + 2.0*x1 + 1.0*x2 + 0.5*x3 + 1.0*x1^2 + 0.5*x2^2 + 0.2*x3^2 + 1.5*x7 + 0.6*x8 + 1.2*x9 + 0.7*x7^2 + 1.2*x8^2 + 0.6*x9^3 + 0.9*x1*x7 + 0.7*x2*x8 + 0.5*x3*x9 + 0.5*x7*x8*x9)
    Zp1 = ex1 / (1 + ex1 + ex2)
    Zp2 = ex2 / (1 + ex1 + ex2)
    Zp3 = 1 - Zp1 - Zp2
    
    Zp1 = ifelse(Zp1 < 10^(-10), 0, Zp1)
    Zp2 = ifelse(Zp2 < 10^(-10), 0, Zp2)
    Zp3 = ifelse(Zp3 < 10^(-10), 0, Zp3)
  }
  Z = NULL
  for (i in 1:n) Z[i] = sample(c(0,1,2), size=1, replace=T, prob=c(Zp1[i],Zp2[i],Zp3[i]))
  
  #table(Z)
  
  if (Aligned==TRUE) {
    # Aligned Y-Model
    if (Ymodel=="Non-Parallel") {
      # Non-parallel
      Yp1 = inv_logit(0.8*x1 + 0.6*x2 + 0.8*x3^2 + 0.5*x7 + 0.7*x8 + 0.5*x9^3 + 0.6*x1*x7)/5
      Yp2 = inv_logit(0.4*x2 + 0.8*x3 + 1.0*x1^2 + 0.3*x8 + 0.3*x9 + 1.2*x7^2 + 0.6*x2*x8)/5
      Yp3 = inv_logit(0.6*x3 + 0.4*x1^2 + 1.0*x2^2 + 0.5*x9 + 0.5*x7^2 + 1.0*x8^2 + 0.6*x3*x9)/5
    } else if (Ymodel=="Parallel") {
      # Parallel
      tau=0.05
      Yp1 = inv_logit(0.8*x1 + 0.6*x2 + 0.8*x3^2 + 0.5*x7 + 0.7*x8 + 0.5*x9^3 + 0.6*x1*x7)/5
      Yp2 = inv_logit(0.8*x1 + 0.6*x2 + 0.8*x3^2 + 0.5*x7 + 0.7*x8 + 0.5*x9^3 + 0.6*x1*x7)/5 + tau
      Yp3 = inv_logit(0.8*x1 + 0.6*x2 + 0.8*x3^2 + 0.5*x7 + 0.7*x8 + 0.5*x9^3 + 0.6*x1*x7)/5 + 2*tau
    }
    
    # potential outcomes
    Y1 = Y2 = Y3= NULL
    for (i in 1:n) {
      Y1[i] = rbinom(1, 1, Yp1[i])
      Y2[i] = rbinom(1, 1, Yp2[i])
      Y3[i] = rbinom(1, 1, Yp3[i])
    }
    
  } else {
    # Not as aligned Y-Model
    x4 = X[,4]; x5 = X[,5]; x6 = X[,6]
    
    if (Ymodel=="Non-Parallel") {
      # Non-parallel
      Yp1 = inv_logit(0.8*x1 + 0.6*x2 + 0.8*x3^2 + 0.8*x4 + 0.7*x5 + 1.2*x6 + 0.5*x4*x6 + 0.5*x9 + 0.2*x9^3)/5
      Yp2 = inv_logit(0.4*x2 + 0.8*x3 + 1.0*x1^2 + 0.6*x5 + 0.5*x6 + 1.0*x4^2 + 0.5*x5*x6 + 1.2*x9^3)/5
      Yp3 = inv_logit(0.6*x3 + 0.4*x1^2 + 1.0*x2^2 + 0.8*x6 + 1.0*x4^2 + 1.2*x5^2 + 1.0*x3*x9)/5
    } else if (Ymodel=="Parallel") {
      # Parallel
      tau=0.05
      Yp1 = inv_logit(0.8*x1 + 0.6*x2 + 0.8*x3^2 + 0.8*x4 + 0.7*x5 + 1.2*x6 + 0.5*x4*x6 + 0.5*x9 + 0.2*x9^3)/5
      Yp2 = inv_logit(0.8*x1 + 0.6*x2 + 0.8*x3^2 + 0.8*x4 + 0.7*x5 + 1.2*x6 + 0.5*x4*x6 + 0.5*x9 + 0.2*x9^3)/5 + tau
      Yp3 = inv_logit(0.8*x1 + 0.6*x2 + 0.8*x3^2 + 0.8*x4 + 0.7*x5 + 1.2*x6 + 0.5*x4*x6 + 0.5*x9 + 0.2*x9^3)/5 + 2*tau
    }
    
    Y1 = Y2 = Y3= NULL
    for (i in 1:n) {
      Y1[i] = rbinom(1, 1, Yp1[i])
      Y2[i] = rbinom(1, 1, Yp2[i])
      Y3[i] = rbinom(1, 1, Yp3[i])
    }
  }
  
  # observed outcomes
  Y = cbind(Y1,Y2,Y3)
  YZ = cbind(Y,Z)
  Yobs = apply(YZ, 1, function(x) x[1:3][x[4]+1])
  
  # ATE(1,2), ATE(1,3), ATE(2,3)
  ATE12 = mean(Y[,1]) - mean(Y[,2])
  ATE13 = mean(Y[,1]) - mean(Y[,3])
  ATE23 = mean(Y[,2]) - mean(Y[,3])
  
  # ATT(1,2), ATT(1,3)
  ATT12 = mean(Y[Z==0,1]) - mean(Y[Z==0,2])
  ATT13 = mean(Y[Z==0,1]) - mean(Y[Z==0,3])
  ATT23 = mean(Y[Z==0,2]) - mean(Y[Z==0,3])
  
  df.sum <- data.frame(n, ratio, p, Zmodel, Ymodel, Aligned, ATE12, ATE13, ATE23, ATT12, ATT13, ATT23)
  data.overall <- data.frame(Y, Yobs, Z, X)
  return(list(data.overall, df.sum))
}



## configuration 0
n <- 100000

## configuration 1
ratio.sim <- c(0.5, 1, 2)

## configuration 2
p.sim <- c(10, 50)

## configuration 3
Zmodel.sim <- c("Linear", "Nonlinear")

## configuration 4
Ymodel.sim <- c("Parallel", "Non-Parallel")

## configuration 5
Aligned.sim <- c(TRUE, FALSE)

set.seed(2018)

df.out <- list()
sim.ticker <- 1
for (i in 1:3){
  for (j in 1:2){
    for(k in 1:2){
      for(l in 1:2){
        for(m in 1:2){
          df.out[[sim.ticker]] <- data_generation_all(n = n, ratio.sim[i], p.sim[j], Zmodel.sim[k], Ymodel.sim[l], Aligned.sim[m])
          print(c(i, j, k, l, m))
          sim.ticker <- sim.ticker + 1
        }
      }
    }
  }
}




save(df.out, file="~/Dropbox/ci-bart/Data/Simulation_complete")


