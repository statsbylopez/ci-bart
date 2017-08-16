rm(list=ls())
source("config.R")

#---------------#
# SET UP DATA   #
#---------------#

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

#proportion of death in each trx group 
group_by(dat,trx) %>% summarise(n=n(), p = sum(death)/n)

#define outcome 5-year survival
#variable: death_5y=1 if died before 5 years since diagnosis
dat<- dat %>% mutate(death_5y = ifelse(death==1 & fupmonth<=60,1,0))


