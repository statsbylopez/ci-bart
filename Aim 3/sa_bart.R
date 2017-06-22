
rm(list=ls())
source("config.R")
library(BayesTree)
igma = 1.0 #y = f(x) + sigma*z , z~N(0,1)
n = 100 #number of observations
set.seed(99)
x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
f= function(x){
  5*(x[,1]*x[,2]) + 2*(x[,3]-.5)^2+2*x[,4]-5*x[,5]
}

Eys = f(x)
prob<-exp(Eys)/(1+exp(Eys))
prob

runis <- runif(n,0,1)
Y <- ifelse(runis < prob,1,0)
bartFit = bart(x,Y,ndpost=200)

plot(bartFit)

fitmat = cbind(prob,bartFit$yhat.train.mean)
colnames(fitmat) = c('y','Ey','lm','bart')
print(cor(fitmat))

#adjust binary offset = sensitivity analysis 
a<-rep(c(qnorm(.2),qnorm(.5),qnorm(.3),qnorm(.4)),each=25)
bartFit2 = bart(x,Y,ndpost=200,binaryOffset =a)
plot(bartFit2)


#--------------------------------------#
# SET UP AND RUN BART 
#------------------------------------#

dat0 <- read.csv(file.path(root, "dat1.csv"))
#some patients are treated at diagnosis
rx0<-filter(dat0,days2trx==0) #548 patients
dat0[dat0$days2trx==0,]$days2trx<-.5
#which variable has missing data
colSums(is.na(dat0))
dat<-na.omit(dat0)


dat[dat$CS_SITESPECIFIC_FACTOR_8<=6,]$CS_SITESPECIFIC_FACTOR_8<-6
dat[!dat$RACE %in% c(1,2),]$RACE<-3 #combine non-white & non-black race group
dat[dat$SPANISH_HISPANIC_ORIGIN %in% c(1:8),]$SPANISH_HISPANIC_ORIGIN<- 1
dat<-dat %>% mutate (tstage =ifelse(TNM_CLIN_T %in% c("3","3A","3B","4"),2,1))%>% select(-TNM_CLIN_T)

names(dat)<-c("age","race","spanish","insurance","income","education",
              "deyo","dyear","psa","gs", "surgsite", "regdose", "boostdose","surgseq","hormone",
              "fupmonth","id","trx","days2trx","death","totdose","tstage")

dat<- dat %>% mutate(dyearcat=ifelse(dyear %in% c(2004:2007), 1, ifelse(dyear %in% c(2008:2010),2,3)))
dat$dyearcat2<- as.factor(dat$dyear)

which(sapply(dat,is.factor))
numindex<-which(sapply(dat,is.numeric)) #age and psa are continous
numindex

ind2<-numindex[c("race","spanish","insurance","income","education",
                 "deyo","gs","trx","tstage","dyearcat")]

facVars = lapply(dat[ , ind2],
                 function(x) {
                   x = as.factor(x)
                   x
                 })

dat1<-cbind(facVars, dat[ , - ind2])

which(sapply(dat1,is.factor))






