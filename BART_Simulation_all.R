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
sample.size <- c(700, 4900)
dat <- lapply(configurations, get_config) %>% bind_rows() 

dat1 <- dat %>% mutate(sample.size = sample.size[1])
dat2 <- dat %>% mutate(sample.size = sample.size[2])

dat <- bind_rows(dat1, dat2) %>% mutate(config.id = 1:n())
#dat.interesting <- dat %>% filter(p == 10, ratio == 2, sample.size == 700)
dat.interesting <- dat

nsim <- 3


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
    scenario.sim <- cbind(scenario.sim, temp)
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
    
    
    ### Vector matching
    n0 <- sum(scenario.sim$Z == 1)
    clustnum <- 5
    
    temp.2 <- kmeans(logit(scenario.sim$p.2), clustnum)
    scenario.sim$Quint.2<-temp.2$cluster
    
    temp.3 <- kmeans(logit(scenario.sim$p.3), clustnum)
    scenario.sim$Quint.3 <- temp.3$cluster
    
    
    ## subgroups to do the matching
    temp.2 <- filter(scenario.sim, Z != 3)  ## where we match r to s
    temp.3 <- filter(scenario.sim, Z != 2)  ## where we match r to t
    
    
    m.2 <- Matchby(Y = temp.2$Yobs, Tr = temp.2$Z == 1, 
                   X = logit(temp.2$p.1), by = temp.2$Quint.3, 
                   calip = 0.5*sd(logit(temp.2$p.1)), replace = T, estimand = "ATT", print.level = 0) 
    
    m.3 <- Matchby(Y = temp.3$Yobs, Tr = temp.3$Z == 1, 
                   X = logit(temp.3$p.1), by = temp.3$Quint.2, 
                   calip = 0.5*sd(logit(temp.3$p.1)), replace = T, estimand = "ATT", print.level = 0) 
    
    ### Identify the matched subgroups 
    rownames(scenario.sim) <- 1:nrow(scenario.sim)
    scenario.sim$id <- 1:nrow(scenario.sim)
    scenario.sim$both <- scenario.sim$id %in% m.2$index.treated & scenario.sim$id %in% m.3$index.treated
    
    m.temp <- scenario.sim[scenario.sim$both == "TRUE", ]
    
    scenario.sim$match.2 <- NULL
    match.12 <- cbind(m.2$index.treated, m.2$index.control)
    temp.2[match.12[1,],]  
    
    match.13 <- cbind(m.3$index.treated, m.3$index.control)
    temp.3[match.13[1,],] 
    
    match.13[,2] <- match.13[,2] + sum(scenario.sim$Z == 2)
    scenario.sim[match.13[10],]
    
    
    match.12 <- match.12[match.12[,1] %in% rownames(m.temp), ]
    match.13 <- match.13[match.13[,1] %in% rownames(m.temp), ]
    
    triplets <- cbind(match.12[order(match.12[,1]), ], match.13[order(match.13[,1]), ])
    triplets <- as.matrix(triplets[,c(1,2,4)])
    
    ## Check for triplets to be well matched on p.0, p.1, p.2
    scenario.sim[triplets[30,],]
    scenario.sim[triplets[31,],]
    
    df.triplets <- rbind(scenario.sim[as.vector(t(triplets)), ])
    df.matched <- rbind(scenario.sim[triplets[,1], ], scenario.sim[triplets[,2], ], scenario.sim[triplets[,3], ])
    percent.matched <- nrow(triplets)/sum(scenario.sim$Z == 1)

    # Matching Estimator
    # For subjects receiving reference treatment (Nom.Treatment = 0.None)
    Yr.hat <- scenario.sim$Yobs[triplets[,1]]
    Ys.hat <- scenario.sim$Yobs[triplets[,2]]
    Yt.hat <- scenario.sim$Yobs[triplets[,3]]

    att12.vm <- mean(Yr.hat - Ys.hat)
    att13.vm <- mean(Yr.hat - Yt.hat)
    #tau23.1.vm <- mean(Ys.hat - Yt.hat)
    
    att12.vm
    m.2$est ### should be near the VM estimate (some people are not matched)
    
    att13.vm
    m.3$est ### should be near the VM estimate (some people are not matched)
    
    
    ### IPTW estimates
    
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
    
    
    ## Regression adjustment (Model-based multiple imputation) 
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
    
    
    ## BART
    
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
    
    
    
    
    
    ## BART_CV (cross-validation, K=5)
    
    xt = scenario.sim[,5:ncol(scenario.sim)]
    
    # Create five equally size folds
    K = 5
    folds = cut(seq(1,nrow(xt)), breaks=K, labels=FALSE)
    
    hyperpar = 1:8 # hyperparameters for k
    cv.error = NULL
    for (s in 1:8) {
      # Perform five-fold cross validation
      mis.error = NULL
      for (k in 1:K) {
        # Segement your data by fold using the which() function
        testIndexes = which(folds==k, arr.ind=TRUE)
        testX = xt[testIndexes, ]
        trainX = xt[-testIndexes, ]
        testY = Y[testIndexes]
        trainY = Y[-testIndexes]
        
        # fit BART
        n.test = length(testY)
        bart_mod = pbart(x.train = trainX, y.train = trainY, k=s, ntree=100, ndpost=n.imps, nskip=500, printevery=2000L)
        bart_pred = pwbart(testX, bart_mod$treedraws)
        y.hat = rbinom(n.test, 1, pnorm(colMeans(bart_pred)))
        
        # calculate misclassfication rate
        mis.error[k] = sum(testY != y.hat) / n.test
      }
      # calculate cross-validation error
      cv.error[s] = mean(mis.error)
    }
    
    # fit the final selected BART model
    best.mod[i] = which.min(cv.error)
    bart_mod = pbart(x.train = xt, y.train = Y, k=best.mod[i], ntree=100, printevery=2000L, ndpost=n.imps, nskip=500)
    
    # treatment 1 to 2 or 3
    xp1 = xt[Z==1,]
    xp2 = xp1
    xp3 = xp1
    xp2[,1] = 2  # switch treatment label 1 to 2
    xp3[,1] = 3  # switch treatment label 1 to 3
    
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
    
    att12.bart.cv = mean(att12.est)
    att13.bart.cv = mean(att13.est)
    
    
    
    
    
    tau.sim <- data.frame(att12.vm, att13.vm, att12.ipw, att13.ipw, att12.ra, att13.ra, att12.bart, att13.bart, att12.bart.cv, att13.bart.cv)
    tau.all <- rbind(tau.all, tau.sim)
    }
    
  tau.all <- tau.all %>% mutate(sim.id = 1:n())
  config.summary <- dat.interesting %>% filter(config.id == config.sim)  
  config.out <- data.frame(config.summary, tau.all)
  config.all <- bind_rows(config.all, config.out)
}

####
write.csv(config.all, "~/SimResults_All.csv", row.names = FALSE)

###
config.all <- read_csv("~/SimResults_All.csv")

config.all.long <- gather(config.all, "type", "estimate", att12.vm:att13.bart) %>% 
  mutate(method = case_when(type == "att12.vm"|type == "att13.vm" ~ "Matching", 
                            type == "att12.ipw"|type == "att13.ipw" ~ "IPTW", 
                            type == "att12.ra"|type == "att13.ra" ~ "RA", 
                            type == "att12.bart"|type == "att13.bart" ~ "BART"), 
         treatment.effect = ifelse(type %in% c("att12.vm", "att12.ipw", "att12.ra", "att12.bart"), "Treatment 2 v. 1", "Treatment 3 v. 1"), 
         true.effect = ifelse(treatment.effect == "Treatment 2 v. 1", ATT12, ATT13), 
         bias = estimate - true.effect, 
         scenario = paste(ratio, p, Zmodel, Ymodel, Aligned, sample.size))

ggplot(config.all.long, aes(x = method, y = bias, fill = method)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = 0), lty = 2, col = "red") +
  scale_fill_brewer(type = "qual") + 
  #scale_color_discrete("Method") + xlim(c(-0.12, 0.12)) + 
  ylab("Bias") + xlab("") + 
  coord_flip() + 
  theme(legend.position="none") + 
  facet_wrap(~scenario + treatment.effect, ncol = 4)

config.all.long.summary <- config.all.long %>% 
  group_by(method, scenario, treatment.effect) %>% 
  summarise(mean.bias = round(mean(bias), 3), 
            rmse = round(sqrt(mean((bias)^2)), 3))
config.all.long.summary

config.all.long.summary %>% 
  ggplot(aes(mean.bias, rmse, colour = method)) + geom_point() + 
  xlab("Average bias") + ylab("Average RMSE") + facet_wrap(~treatment.effect)


########################################################################
###### ANOVA's
########################################################################
library(broom)
anova.matching <- lm(bias~(Zmodel + Ymodel + Aligned + type)^2-1, data = filter(config.all.long, method == "Matching"))
anova.fit.matching <- tidy(anova(anova.matching)) %>% arrange(-meansq) %>% 
  mutate(fraction.explained = meansq/sum(meansq), 
         type = "Matching", rank = 1:n())

anova.iptw <- lm(bias~(Zmodel + Ymodel + Aligned + type)^2-1, data = filter(config.all.long, method == "IPTW"))
anova.fit.iptw <- tidy(anova(anova.iptw)) %>% arrange(-meansq) %>% 
  mutate(fraction.explained = meansq/sum(meansq), 
         type = "IPTW", rank = 1:n())

anova.bart <- lm(bias~(Zmodel + Ymodel + Aligned + type)^2-1, data = filter(config.all.long, method == "BART"))
anova.fit.bart <- tidy(anova(anova.bart)) %>% arrange(-meansq) %>% 
  mutate(fraction.explained = meansq/sum(meansq), 
         type = "BART", rank = 1:n())

anova.ra <- lm(bias~(Zmodel + Ymodel + Aligned + type)^2-1, data = filter(config.all.long, method == "RA"))
anova.fit.ra <- tidy(anova(anova.ra)) %>% arrange(-meansq) %>% 
  mutate(fraction.explained = meansq/sum(meansq), 
         type = "RA", rank = 1:n())


bind_rows(anova.fit.matching, anova.fit.iptw, anova.fit.bart, anova.fit.ra) %>% 
  filter(rank < 6) %>%
  ggplot(aes(fraction.explained, x = type, fill = reorder(term, fraction.explained))) + 
  geom_col() + 
  scale_fill_brewer(type = "qual", "Variable") + 
  labs(title = "Percent of bias explained by each configuration factor", subtitle = "Top factors only")



########################################################################
###### Extras from last time
########################################################################

load("Data/Results_BART_Scenario1.RData")
scen1.bart <- data.frame(att.est)
scen1.bart$scenario <- "Aligned"
scen1.bart <- scen1.bart %>% gather("type", "estimate", att12:att13) %>% 
  mutate(method = "BART", 
         treatment.effect = ifelse(type == "att12", "Treatment 2 v. 1", "Treatment 3 v 1")) 

load("Data/Results_RA_Scenario1.RData")
scen1.ra <- data.frame(att.est)
scen1.ra$scenario <- "Aligned"
scen1.ra <- scen1.ra %>% gather("type", "estimate", att12:att13) %>% 
  mutate(method = "RA", 
         treatment.effect = ifelse(type == "att12", "Treatment 2 v. 1", "Treatment 3 v 1")) 



load("Data/Results_BART_Scenario2.RData")
scen2.bart <- data.frame(att.est)
scen2.bart$scenario <- "Less Aligned"
scen2.bart <- scen2.bart %>% gather("type", "estimate", att12:att13) %>% 
  mutate(method = "BART", 
         treatment.effect = ifelse(type == "att12", "Treatment 2 v. 1", "Treatment 3 v 1")) 

load("Data/Results_RA_Scenario2.RData")
scen2.ra <- data.frame(att.est)
scen2.ra$scenario <- "Less Aligned"
scen2.ra <- scen2.ra %>% gather("type", "estimate", att12:att13) %>% 
  mutate(method = "RA",  
         treatment.effect = ifelse(type == "att12", "Treatment 2 v. 1", "Treatment 3 v 1")) 

data.out.long.all <- bind_rows(data.out.long, 
                               scen1.bart, scen1.ra, scen2.bart, scen2.ra)

data.out.long.all <- data.out.long.all %>% left_join(true.ate)

ggplot(data.out.long.all, aes(estimate, colour = method)) + 
  geom_density() + 
  geom_vline(data = true.ate, aes(xintercept = effect), lty = 2, col = "black") +
  scale_color_discrete("Method") + xlim(c(-0.2, 0.05)) + 
  facet_wrap(~scenario + treatment.effect)


ggplot(data.out.long.all, aes(estimate - effect, colour = method)) + 
  geom_density() + 
  geom_vline(aes(xintercept = 0), lty = 2, col = "black") +
  scale_color_discrete("Method") + xlim(c(-0.12, 0.12)) + 
  xlab("Bias") + ylab("Density") + 
  facet_wrap(~scenario + treatment.effect)

ggplot(data.out.long.all, aes(x = method, y = estimate - effect, fill = method)) + 
  geom_boxplot() + 
  #geom_vline(aes(xintercept = 0), lty = 2, col = "black") +
  #scale_color_discrete("Method") + xlim(c(-0.12, 0.12)) + 
  ylab("Bias") + xlab("") + 
  coord_flip() + ylim(-0.2, 0.2) + 
  theme(legend.position="none") + 
  facet_wrap(~scenario + treatment.effect, scales = "free")

data.out.long.all %>% group_by(method, scenario, treatment.effect) %>% 
  summarise(mean.bias = round(mean(estimate - effect), 3), 
            rmse = round(sqrt(mean((estimate - effect)^2)), 3))

