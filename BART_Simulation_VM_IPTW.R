library(nnet) 
library(tidyverse)
library(Matching)


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
dat.interesting <- dat %>% filter(p == 10, ratio == 2, sample.size == 700)

nsim <- 5


####################################################################
####################################################################
#### Vector Matching & IPTW
####################################################################
####################################################################


### Everyone 
set.seed(0)
config.all <- NULL ### Stores all of the simulations

for (j in 1:nrow(dat.interesting)){
  config.sim <- dat.interesting$config.id[j]
  scenario <- data.frame(df.out[[config.sim]][1]) %>% mutate(Z = Z + 1)
  tau.all <- NULL
  
    for (i in 1:nsim){
      
    print(paste("simulation", i, "with scenario", j))
    
    scenario.sim <- sample_n(scenario, 3500)
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
                   calip = 0.5*sd(logit(temp.2$p.1)), replace = T, estimand = "ATT") 
    
    m.3 <- Matchby(Y = temp.3$Yobs, Tr = temp.3$Z == 1, 
                   X = logit(temp.3$p.1), by = temp.3$Quint.2, 
                   calip = 0.5*sd(logit(temp.3$p.1)), replace = T, estimand = "ATT") 
    
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
    scenario.sim[match.13[100,],]
    
    
    match.12 <- match.12[match.12[,1] %in% rownames(m.temp), ]
    match.13 <- match.13[match.13[,1] %in% rownames(m.temp), ]
    
    triplets <- cbind(match.12[order(match.12[,1]), ], match.13[order(match.13[,1]), ])
    triplets <- as.matrix(triplets[,c(1,2,4)])
    
    ## Check for triplets to be well matched on p.0, p.1, p.2
    scenario.sim[triplets[30,],]
    scenario.sim[triplets[31,],]
    
    df.triplets <- rbind(scenario.sim[as.vector(t(triplets)), ])
    df.matched <- rbind(scenario.sim[triplets[,1], ], scenario.sim[triplets[,2], ], scenario.sim[triplets[,3], ])
    
    df.unique <- unique(df.matched)
    
    percent.matched <- nrow(triplets)/sum(scenario.sim$Z == 1)
    percent.matched
    
    
    # Matching Estimator
    # For subjects receiving reference treatment (Nom.Treatment = 0.None)
    Yr.hat <- scenario.sim$Yobs[triplets[,1]]
    Ys.hat <- scenario.sim$Yobs[triplets[,2]]
    Yt.hat <- scenario.sim$Yobs[triplets[,3]]
    
    #scenario.sim$Tr <- scenario.sim$Z == 0
    #scenario.sim$Ts <- scenario.sim$Z == 1
    #scenario.sim$Tt <- scenario.sim$Z == 2
    #
    #n0m <- sum(scenario.sim$Tr)
    #n1m <- sum(scenario.sim$Ts)
    #n2m <- sum(scenario.sim$Tt)
    #Nm <- n0m + n1m + n2m
    #
    #
    ### Number of times matched
    #triplets <- data.frame(id = as.vector(triplets))
    #triplets.count <- triplets %>% group_by(id) %>% count() %>% rename(psi1 = n)
    #triplets.count %>% arrange(-psi1) %>% head()
    #scenario.sim <- scenario.sim %>% left_join(triplets.count)
    #scenario.sim[is.na(scenario.sim$psi1),]$psi1 <- 0
    
    
    tau12.1.vm <- mean(Yr.hat - Ys.hat)
    tau13.1.vm <- mean(Yr.hat - Yt.hat)
    tau23.1.vm <- mean(Ys.hat - Yt.hat)
    
    tau12.1.vm
    m.2$est ### should be near the VM estimate (some people are not matched)
    
    tau13.1.vm
    m.3$est ### should be near the VM estimate (some people are not matched)
    
    
    ### IPTW estimates
    
    mu1.hat.ratio.logit = sum(y[treat==1] / pred.class.probs.logit[treat==1,1]) / sum(1 / pred.class.probs.logit[treat==1,1])
    mu2.hat.ratio.logit = sum(y[treat==2] / pred.class.probs.logit[treat==2,2]) / sum(1 / pred.class.probs.logit[treat==2,2])
    mu3.hat.ratio.logit = sum(y[treat==3] / pred.class.probs.logit[treat==3,3]) / sum(1 / pred.class.probs.logit[treat==3,3])
    
    tau12.1.ipw = mu1.hat.ratio.logit - mu2.hat.ratio.logit
    tau13.1.ipw = mu1.hat.ratio.logit - mu3.hat.ratio.logit
    tau23.1.ipw = mu2.hat.ratio.logit - mu3.hat.ratio.logit
    
    
    tau.sim <- data.frame(tau12.1.vm, tau13.1.vm, tau23.1.vm, tau12.1.ipw, tau13.1.ipw, tau23.1.ipw)
    tau.all <- rbind(tau.all, tau.sim)
    }
    
  tau.all <- tau.all %>% mutate(sim.id = 1:n())
  config.summary <- dat.interesting %>% filter(config.id == config.sim)  
  config.out <- data.frame(config.summary, tau.all)
  config.all <- bind_rows(config.all, config.out)
}

  
#data.out <- bind_rows(scen1.out, scen2.out)
#true.ate <- data.frame(treatment.effect = c("Treatment 2 v. 1", "Treatment 3 v 1", "Treatment 2 v. 1", "Treatment 3 v 1"), 
#                       scenario = c("Aligned", "Aligned", "Less Aligned", "Less Aligned"), 
#                       effect = c(mydata1$Ests$ATT12, mydata1$Ests$ATT13, mydata2$Ests$ATT12, mydata2$Ests$ATT13))

config.all.long <- gather(config.all, "type", "estimate", tau12.1.vm:tau23.1.ipw) %>% 
  mutate(method = ifelse(type == "tau12.1.vm" | type == "tau13.1.vm"| type == "tau23.1.vm", "Matching", "IPTW"), 
         treatment.effect = ifelse(type == "tau12.1.vm" | type == "tau12.1.ipw", "Treatment 2 v. 1", 
                                   ifelse(type == "tau23.1.vm"|type == "tau23.1.ipw", "Treatment 3 v. 2", "Treatment 3 v. 1")), 
         true.effect = ifelse(treatment.effect == "Treatment 2 v. 1", ATT12, 
                              ifelse(treatment.effect == "Treatment 3 v. 1", ATT13, ATT23)), 
         bias = estimate - true.effect, 
         scenario = paste(ratio, p, Zmodel, Ymodel, Aligned))

ggplot(config.all.long, aes(x = method, y = bias, fill = method)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = 0), lty = 2, col = "red") +
  scale_fill_brewer(type = "qual") + 
  #scale_color_discrete("Method") + xlim(c(-0.12, 0.12)) + 
  ylab("Bias") + xlab("") + 
  coord_flip() + 
  theme(legend.position="none") + 
  facet_wrap(~scenario + treatment.effect)


config.all.long %>% group_by(method, scenario, treatment.effect) %>% 
  summarise(mean.bias = round(mean(bias), 3), 
            rmse = round(sqrt(mean((bias)^2)), 3))



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

