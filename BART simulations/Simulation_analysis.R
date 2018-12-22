library(nnet) 
library(tidyverse)
library(Matching)
library(arm)
#library(dbarts)
library(BART)

####################################################################
####################################################################
#### Estimates from simulation
####################################################################
####################################################################
set.seed(0)
logit <- function(p){return(log(p/(1-p)))}

load(file="~/Dropbox/ci-bart/Data/Simulation_complete")

configurations <- 1:48
get_config <- function(config) {
  message(paste("reading", config, "data..."))
  summary.config <- data.frame(df.out[[config]][2]) 
  return(summary.config)
}
get_config(25)


# Extract the values we will need
sample.size <- c(700, 3500)
dat <- lapply(configurations, get_config) %>% bind_rows() 

dat1 <- dat %>% mutate(sample.size = sample.size[1])
dat2 <- dat %>% mutate(sample.size = sample.size[2])

dat <- bind_rows(dat1, dat2) %>% mutate(config.id = 1:n())
#dat.interesting <- dat %>% filter(p == 10, ratio == 2, sample.size == 700)
dat.interesting <- dat

nsim <- 200
ptm <- proc.time()

####################################################################
####################################################################
#### Vector Matching & IPTW
####################################################################
####################################################################

config.sim <- 96
ifelse(config.sim <=48, config.sim, config.sim - 48)

### Everyone 
set.seed(0)
config.all <- NULL ### Stores all of the simulations

for (j in 1:nrow(dat.interesting)){
  config.sim <- dat.interesting$config.id[j]
  sample.size <- dat.interesting$sample.size[j]
  config.setting <- ifelse(config.sim <=48, config.sim, config.sim - 48)
  scenario <- data.frame(df.out[[config.setting]][1]) %>% mutate(Z = Z + 1)
  tau.all <- NULL
  best.mod = NULL
  
    for (i in 1:nsim){
      
    print(paste("simulation", i, "with scenario", j))
    
    scenario.sim <- sample_n(scenario, sample.size)  ### Make sure this works with 4900 
    scenario.sim <- arrange(scenario.sim, Z)
    y <- scenario.sim$Yobs
    treat <- scenario.sim$Z
    mlogit.data <- scenario.sim[,5:ncol(scenario.sim)]
    
    fit1 <- multinom(Z ~ ., data = mlogit.data, trace = FALSE)
    pred.class.probs.logit <- fitted(fit1)
    temp <- data.frame(fitted(fit1))
    colnames(temp) <- c("p.1", "p.2", "p.3")
    scenario.sim.vm <- cbind(scenario.sim, temp)
    #
    #
    #### Drop subjects with extreme GPS & refit model
    #min.max.Ps <- scenario.sim %>%
    #  group_by(Z) %>%
    #  summarize(min0 = min(p.0), max0 = max(p.0), 
    #            min1 = min(p.1), max1 = max(p.1), 
    #            min2 = min(p.2), max2 = max(p.2))
    #
    #scenario.sim$Eligible <- 
    #  scenario.sim$p.0 >= max(min.max.Ps$min0) & scenario.sim$p.0 <= min(min.max.Ps$max0) &
    #  scenario.sim$p.1 >= max(min.max.Ps$min1) & scenario.sim$p.1 <= min(min.max.Ps$max1) &
    #  scenario.sim$p.2 >= max(min.max.Ps$min2) & scenario.sim$p.2 <= min(min.max.Ps$max2)
    #scenario.sim <- filter(scenario.sim, Eligible)
    #
    #fit2 <- multinom(Z ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7 + X.8 + X.9 + X.10, data = scenario.sim, trace = FALSE)
    #temp <- data.frame(fitted(fit2))
    #colnames(temp) <- c("p.0", "p.1", "p.2")
    #scenario.sim <- scenario.sim %>%
    #  dplyr::select(-p.0, -p.1, -p.2)
    #scenario.sim <- cbind(scenario.sim, temp)
    
    
    ###################################
    ## Matching on a vector
    ###################################
    
    n0 <- sum(scenario.sim.vm$Z == 1)
    clustnum <- 5
    
    temp.2 <- kmeans(logit(scenario.sim.vm$p.2), clustnum)
    scenario.sim.vm$Quint.2<-temp.2$cluster
    
    temp.3 <- kmeans(logit(scenario.sim.vm$p.3), clustnum)
    scenario.sim.vm$Quint.3 <- temp.3$cluster
    
    
    ## subgroups to do the matching
    temp.2 <- filter(scenario.sim.vm, Z != 3)  ## where we match r to s
    temp.3 <- filter(scenario.sim.vm, Z != 2)  ## where we match r to t
    
    
    m.2 <- Matchby(Y = temp.2$Yobs, Tr = temp.2$Z == 1, 
                   X = logit(temp.2$p.1), by = temp.2$Quint.3, 
                   calip = 0.5*sd(logit(temp.2$p.1)), replace = T, estimand = "ATT", print.level = 0) 
    
    m.3 <- Matchby(Y = temp.3$Yobs, Tr = temp.3$Z == 1, 
                   X = logit(temp.3$p.1), by = temp.3$Quint.2, 
                   calip = 0.5*sd(logit(temp.3$p.1)), replace = T, estimand = "ATT", print.level = 0) 
    
    ### Identify the matched subgroups 
    rownames(scenario.sim.vm) <- 1:nrow(scenario.sim.vm)
    scenario.sim.vm$id <- 1:nrow(scenario.sim.vm)
    scenario.sim.vm$both <- scenario.sim.vm$id %in% m.2$index.treated & scenario.sim.vm$id %in% m.3$index.treated
    
    m.temp <- scenario.sim.vm[scenario.sim.vm$both == "TRUE", ]
    
    scenario.sim.vm$match.2 <- NULL
    match.12 <- cbind(m.2$index.treated, m.2$index.control)
    temp.2[match.12[1,],]  
    
    match.13 <- cbind(m.3$index.treated, m.3$index.control)
    temp.3[match.13[1,],] 
    
    match.13[,2] <- match.13[,2] + sum(scenario.sim.vm$Z == 2)
    scenario.sim[match.13[10],]
    
    
    match.12 <- match.12[match.12[,1] %in% rownames(m.temp), ]
    match.13 <- match.13[match.13[,1] %in% rownames(m.temp), ]
    
    triplets <- cbind(match.12[order(match.12[,1]), ], match.13[order(match.13[,1]), ])
    triplets <- as.matrix(triplets[,c(1,2,4)])
    
    ## Check for triplets to be well matched on p.0, p.1, p.2
    scenario.sim.vm[triplets[30,],]
    scenario.sim.vm[triplets[31,],]
    
    df.triplets <- rbind(scenario.sim.vm[as.vector(t(triplets)), ])
    df.matched <- rbind(scenario.sim.vm[triplets[,1], ], scenario.sim.vm[triplets[,2], ], scenario.sim.vm[triplets[,3], ])
    percent.matched <- nrow(triplets)/sum(scenario.sim.vm$Z == 1)
    
    # Matching Estimator
    # For subjects receiving reference treatment (Nom.Treatment = 0.None)
    Yr.hat <- scenario.sim.vm$Yobs[triplets[,1]]
    Ys.hat <- scenario.sim.vm$Yobs[triplets[,2]]
    Yt.hat <- scenario.sim.vm$Yobs[triplets[,3]]
    
    att12.vm <- mean(Yr.hat - Ys.hat)
    att13.vm <- mean(Yr.hat - Yt.hat)
    #tau23.1.vm <- mean(Ys.hat - Yt.hat)
    
    att12.vm
    m.2$est ### should be near the VM estimate (some people are not matched)
    
    att13.vm
    m.3$est ### should be near the VM estimate (some people are not matched)
    
    ###################################
    ## iptw
    ###################################
    
    mu_1.1_hat.iptw = mean(y[treat==1])  ## mu_{t',t'}, the unweighted mean of observed outcomes from units assigned reference
    
    wt_1.2_hat = pred.class.probs.logit[,1]/pred.class.probs.logit[,2]
    mu_1.2_hat.iptw = sum(y[treat==2] * wt_1.2_hat[treat==2]) / sum(wt_1.2_hat[treat==2])
      
    
    wt_1.3_hat = pred.class.probs.logit[,1]/pred.class.probs.logit[,3]
    mu_1.3_hat.iptw = sum(y[treat==3] * wt_1.3_hat[treat == 3]) / sum(wt_1.3_hat[treat == 3])
      
    #mu2.hat.ratio.logit = sum(y[treat==2] / pred.class.probs.logit[treat==2,2]) / sum(1 / pred.class.probs.logit[treat==2,2])
    #mu3.hat.ratio.logit = sum(y[treat==3] / pred.class.probs.logit[treat==3,3]) / sum(1 / pred.class.probs.logit[treat==3,3])
    
    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3710547/
    
    att12.ipw = mu_1.1_hat.iptw - mu_1.2_hat.iptw
    att13.ipw = mu_1.1_hat.iptw - mu_1.3_hat.iptw
    #tau23.1.ipw = mu_1.2_hat.iptw - mu_1.3_hat.iptw
    
    
    ###################################
    ## Regression adjustment (model based imputation)
    ###################################
    
    Y <- scenario.sim$Yobs
    Z <- scenario.sim$Z
    data1 <- scenario.sim[,4:ncol(scenario.sim)] %>% filter(Z == 1)
    data2 <- scenario.sim[,4:ncol(scenario.sim)] %>% filter(Z == 2)
    data3 <- scenario.sim[,4:ncol(scenario.sim)] %>% filter(Z == 3)
    
    # outcome model for each treatment level
    # bayesian logistic regression model, default Cauchy prior with scale 2.5
    mod1 = bayesglm(Yobs ~ ., data = data1, family = binomial(link="logit"), x = TRUE)
    mod2 = bayesglm(Yobs ~ ., data = data2, family = binomial(link="logit"), x = TRUE)
    mod3 = bayesglm(Yobs ~ ., data = data3, family = binomial(link="logit"), x = TRUE)
    
    # simulate the uncertainty in the estiamted coefficients
    n.imps = 1000 # number of imputations
    sim1 = sim(mod1, n.sims = n.imps)
    sim2 = sim(mod2, n.sims = n.imps)
    sim3 = sim(mod3, n.sims = n.imps)
    
    sim1.beta = coef(sim1)
    sim2.beta = coef(sim2)
    sim3.beta = coef(sim3)
    
    # predictive simulation using the binomial distribution
    X1.tilde = model.matrix(mod1)
    n1.tilde = nrow(X1.tilde)
    
    # predict potential outcomes
    y11.tilde = array(NA, c(n.imps, n1.tilde))
    y12.tilde = array(NA, c(n.imps, n1.tilde))
    y13.tilde = array(NA, c(n.imps, n1.tilde))
    
    for (s in 1:n.imps) {
      # predict potential outcome Y1 using X1
      p11.tilde = invlogit(X1.tilde %*% sim1.beta[s,])
      y11.tilde[s,] = rbinom(n1.tilde, 1, p11.tilde)
      
      # predict potential outcome Y2 using X1
      p12.tilde = invlogit(X1.tilde %*% sim2.beta[s,])
      y12.tilde[s,] = rbinom(n1.tilde, 1, p12.tilde)
      
      # predict potential outcome Y3 using X1
      p13.tilde = invlogit(X1.tilde %*% sim3.beta[s,])
      y13.tilde[s,] = rbinom(n1.tilde, 1, p13.tilde)
    }
    
    # Average treatment effects on the treatment group 1 (ATTs)
    att12.est = att13.est = NULL
    for (m in 1:n.imps) {
      # potential outcomes for treatment group 1
      y11.hat = mean(y11.tilde[m,])
      y12.hat = mean(y12.tilde[m,])
      y13.hat = mean(y13.tilde[m,])
      
      # risk difference
      att12.est[m] = y11.hat - y12.hat
      att13.est[m] = y11.hat - y13.hat
    }
    att12.ra = mean(att12.est)
    att13.ra = mean(att13.est)
    
    ###################################
    ## traditional BART
    ###################################
    
    # treatment 1 to 2 or 3
    xt = scenario.sim[,5:ncol(scenario.sim)]
    xp1 = xt[Z==1,]
    xp2 = xp1
    xp3 = xp1
    xp2[,1] = 2  # switch treatment label 1 to 2
    xp3[,1] = 3  # switch treatment label 1 to 3
    
    #bart_tot12 = bart(x.train = xt, y.train = Y,  x.test = xp2, ntree = 100, ndpost = n.imps)
    #bart_tot13 = bart(x.train = xt, y.train = Y,  x.test = xp3, ntree = 100, ndpost = n.imps)
    
    bart_mod = pbart(x.train = xt, y.train = Y, k=2, ntree=100, ndpost=n.imps, nskip=500, printevery=2000L)
    bart_pred1 = pwbart(xp1, bart_mod$treedraws)
    bart_pred2 = pwbart(xp2, bart_mod$treedraws)
    bart_pred3 = pwbart(xp3, bart_mod$treedraws)
    
    # Average treatment effects on the treatment group 1 (ATTs)
    n1 = nrow(xp1)
    att12.est = att13.est = NULL
    for (m in 1:n.imps) {
      # potential outcomes for treatment group 1
      y11.hat = rbinom(n1, 1, pnorm(bart_pred1[m,]))
      y12.hat = rbinom(n1, 1, pnorm(bart_pred2[m,]))
      y13.hat = rbinom(n1, 1, pnorm(bart_pred3[m,]))
      
      # risk difference
      att12.est[m] = mean(y11.hat) - mean(y12.hat)
      att13.est[m] = mean(y11.hat) - mean(y13.hat)
    }
    
    att12.bart = mean(att12.est)
    att13.bart = mean(att13.est)
    
    ###################################
    ## psBART
    ###################################
    
   # # estimate propensity scores using multinomial BART (mbart)
   # ps.mod = mbart(xt, Z, ndpost=n.imps, nskip=500, printevery=2000L)
   # idx = seq(1, 3*length(Z), by=3)
   # ps1 = ps.mod$prob.train.mean[idx]
   # ps2 = ps.mod$prob.train.mean[idx+1]
   # ps1.linear = logit(ps1)
   # ps2.linear = logit(ps2)
   # 
   # # incoporate PS into the BART
   # xx = cbind(xt, ps1.linear, ps2.linear)
   # bart_mod = pbart(x.train = xx, y.train = Y, k=2, ntree=100, ndpost=n.imps, nskip=500, printevery=2000L)
   # 
   # # treatment 1 to 2 or 3
   # xp1 = xx[Z==1,]
   # xp2 = xp1
   # xp3 = xp1
   # xp2[,1] = 2  # switch treatment label 1 to 2
   # xp3[,1] = 3  # switch treatment label 1 to 3
   # 
   # bart_pred1 = pwbart(xp1, bart_mod$treedraws)
   # bart_pred2 = pwbart(xp2, bart_mod$treedraws)
   # bart_pred3 = pwbart(xp3, bart_mod$treedraws)
   # 
   # # Average treatment effects on the treatment group 1 (ATTs)
   # n1 = nrow(xp1)
   # att12.est = att13.est = NULL
   # for (m in 1:n.imps) {
   #   # potential outcomes for treatment group 1
   #   y11.hat = rbinom(n1, 1, pnorm(bart_pred1[m,]))
   #   y12.hat = rbinom(n1, 1, pnorm(bart_pred2[m,]))
   #   y13.hat = rbinom(n1, 1, pnorm(bart_pred3[m,]))
   #   
   #   # risk difference
   #   att12.est[m] = mean(y11.hat) - mean(y12.hat)
   #   att13.est[m] = mean(y11.hat) - mean(y13.hat)
   # }
   # 
   # att12.psbart = mean(att12.est)
   # att13.psbart = mean(att13.est)
   # 
   # 
    
    
    
    ## BART_CV (cross-validation, K=5)
    
   # xt = scenario.sim[,5:ncol(scenario.sim)]
   # 
   # # Create five equally size folds
   # K = 5
   # folds = cut(seq(1,nrow(xt)), breaks=K, labels=FALSE)
   # 
   # hyperpar = 1:8 # hyperparameters for k
   # cv.error = NULL
   # for (s in 1:8) {
   #   # Perform five-fold cross validation
   #   mis.error = NULL
   #   for (k in 1:K) {
   #     # Segement your data by fold using the which() function
   #     testIndexes = which(folds==k, arr.ind=TRUE)
   #     testX = xt[testIndexes, ]
   #     trainX = xt[-testIndexes, ]
   #     testY = Y[testIndexes]
   #     trainY = Y[-testIndexes]
   #     
   #     # fit BART
   #     n.test = length(testY)
   #     bart_mod = pbart(x.train = trainX, y.train = trainY, k=s, ntree=100, ndpost=n.imps, nskip=500, printevery=2000L)
   #     bart_pred = pwbart(testX, bart_mod$treedraws)
   #     y.hat = rbinom(n.test, 1, pnorm(colMeans(bart_pred)))
   #     
   #     # calculate misclassfication rate
   #     mis.error[k] = sum(testY != y.hat) / n.test
   #   }
   #   # calculate cross-validation error
   #   cv.error[s] = mean(mis.error)
   # }
   # 
   # # fit the final selected BART model
   # best.mod[i] = which.min(cv.error)
   # bart_mod = pbart(x.train = xt, y.train = Y, k=best.mod[i], ntree=100, printevery=2000L, ndpost=n.imps, nskip=500)
   # 
   # # treatment 1 to 2 or 3
   # xp1 = xt[Z==1,]
   # xp2 = xp1
   # xp3 = xp1
   # xp2[,1] = 2  # switch treatment label 1 to 2
   # xp3[,1] = 3  # switch treatment label 1 to 3
   # 
   # bart_pred1 = pwbart(xp1, bart_mod$treedraws)
   # bart_pred2 = pwbart(xp2, bart_mod$treedraws)
   # bart_pred3 = pwbart(xp3, bart_mod$treedraws)
   # 
   # # Average treatment effects on the treatment group 1 (ATTs)
   # n1 = nrow(xp1)
   # att12.est = att13.est = NULL
   # for (m in 1:n.imps) {
   #   # potential outcomes for treatment group 1
   #   y11.hat = rbinom(n1, 1, pnorm(bart_pred1[m,]))
   #   y12.hat = rbinom(n1, 1, pnorm(bart_pred2[m,]))
   #   y13.hat = rbinom(n1, 1, pnorm(bart_pred3[m,]))
   #   
   #   # risk difference
   #   att12.est[m] = mean(y11.hat) - mean(y12.hat)
   #   att13.est[m] = mean(y11.hat) - mean(y13.hat)
   # }
   # 
   # att12.bart.cv = mean(att12.est)
   # att13.bart.cv = mean(att13.est)
   # 
   # 
    
    
    
    #tau.sim <- data.frame(att12.vm, att13.vm, att12.ipw, att13.ipw, att12.ra, att13.ra, att12.bart, att13.bart, att12.bart.cv, att13.bart.cv)
    #tau.sim <- data.frame(att12.vm, att13.vm, att12.ipw, att13.ipw, att12.ra, att13.ra, att12.bart, att13.bart, att12.psbart, att13.psbart)
    tau.sim <- data.frame(att12.vm, att13.vm, att12.ipw, att13.ipw, att12.ra, att13.ra, att12.bart, att13.bart)
    tau.all <- rbind(tau.all, tau.sim)
    }
    
  tau.all <- tau.all %>% mutate(sim.id = 1:n())
  config.summary <- dat.interesting %>% filter(config.id == config.sim)  
  config.out <- data.frame(config.summary, tau.all)
  config.all <- bind_rows(config.all, config.out)
  file.temp <- paste0("/Users/mlopez1/Dropbox/ci-bart/Simulations/SimResults_All_configuration", j, ".csv")
  write.csv(config.out, file = file.temp, row.names = FALSE)
}

nsim
proc.time() - ptm

####
write.csv(config.all, "~/SimResults_All.csv", row.names = FALSE)
