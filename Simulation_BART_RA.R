library(tidyverse)
library(arm)
library(dbarts)


####################################################################
####################################################################
#### Estimates from simulation
####################################################################
####################################################################
set.seed(2018)
logit <- function(p){return(log(p/(1-p)))}

load(file="Simulation_complete")

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

nsim <- 200


####################################################################
####################################################################
#### BART & RA
####################################################################
####################################################################

### Everyone
set.seed(2018)
config.all <- NULL ### Stores all of the simulations

n.imps = 1000 # number of imputations

for (j in 1:nrow(dat.interesting)){
    config.sim <- dat.interesting$config.id[j]
    sample.size <- dat.interesting$sample.size[j]
    scenario <- data.frame(df.out[[config.sim]][1]) %>% mutate(Z = Z + 1)
    tau.all <- NULL
    
    for (i in 1:nsim){
        
        print(paste("simulation", i, "with scenario", j))
        
        scenario.sim <- sample_n(scenario, sample.size)  ### Make sure this works with 4900
        scenario.sim <- arrange(scenario.sim, Z)
        Y <- scenario.sim$Yobs
        Z <- scenario.sim$Z
        
        
        ## Regression adjustment (Model-based multiple imputation)
        data1 <- scenario.sim[Z==1, 4:ncol(scenario.sim)]
        data2 <- scenario.sim[Z==2, 4:ncol(scenario.sim)]
        data3 <- scenario.sim[Z==3, 4:ncol(scenario.sim)]
        
        # outcome model for each treatment level
        # bayesian logistic regression model, default Cauchy prior with scale 2.5
        mod1 = bayesglm(Yobs ~ ., data = data1, family = binomial(link="logit"), x = TRUE)
        mod2 = bayesglm(Yobs ~ ., data = data2, family = binomial(link="logit"), x = TRUE)
        mod3 = bayesglm(Yobs ~ ., data = data3, family = binomial(link="logit"), x = TRUE)
        
        # simulate the uncertainty in the estiamted coefficients
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
        
        bart_tot12 = bart(x.train = xt, y.train = Y,  x.test = xp2, ntree = 100, ndpost = n.imps)
        bart_tot13 = bart(x.train = xt, y.train = Y,  x.test = xp3, ntree = 100, ndpost = n.imps)
        
        # Average treatment effects on the treatment group 1 (ATTs)
        n1 = nrow(xp1)
        att12.est = att13.est = NULL
        for (m in 1:n.imps) {
            # potential outcomes for treatment group 1
            y11.hat = rbinom(n1, 1, pnorm(bart_tot12$yhat.train[m,Z==1]))
            y12.hat = rbinom(n1, 1, pnorm(bart_tot12$yhat.test[m,]))
            y13.hat = rbinom(n1, 1, pnorm(bart_tot13$yhat.test[m,]))
            
            # risk difference
            att12.est[m] = mean(y11.hat) - mean(y12.hat)
            att13.est[m] = mean(y11.hat) - mean(y13.hat)
        }
        
        att12.bart = mean(att12.est)
        att13.bart = mean(att13.est)
        
        
        tau.sim <- data.frame(att12.ra, att13.ra, att12.bart, att13.bart)
        tau.all <- rbind(tau.all, tau.sim)
        
        cat(i,"\n")
    }
    
    tau.all <- tau.all %>% mutate(sim.id = 1:n())
    config.summary <- dat.interesting %>% filter(config.id == config.sim)
    config.out <- data.frame(config.summary, tau.all)
    config.all <- bind_rows(config.all, config.out)
}


config.all.long <- gather(config.all, "type", "estimate", att12.ra:att13.bart) %>% mutate(method = ifelse(type == "att12.ra" | type == "att13.ra", "RA", "BART"), treatment.effect = ifelse(type == "att12.bart" | type == "att12.ra", "Treatment 1 v. 2", "Treatment 1 v. 3"), true.effect = ifelse(treatment.effect == "Treatment 1 v. 2", ATT12, ifelse(treatment.effect == "Treatment 1 v. 3", ATT13, ATT23)), bias = estimate - true.effect, scenario = paste(ratio, p, Zmodel, Ymodel, Aligned))

ggplot(config.all.long, aes(x = method, y = bias, fill = method)) +
geom_boxplot() +
geom_hline(aes(yintercept = 0), lty = 2, col = "red") +
scale_fill_brewer(type = "qual") +
#scale_color_discrete("Method") + xlim(c(-0.12, 0.12)) +
ylab("Bias") + xlab("") +
coord_flip() +
theme(legend.position="none") +
facet_wrap(~scenario + treatment.effect, ncol = 4)


config.all.long %>% group_by(method, scenario, treatment.effect) %>% summarise(mean.bias = round(mean(bias), 3), rmse = round(sqrt(mean((bias)^2)), 3))






