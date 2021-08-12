## this script fits methods (ours and others) to 25 iterations of simulated data

job=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
a=(25*(job-1))+1
b=a+24

library(dplyr);library(survival);library(splitstackshape);library(splines);library(zoo)

load("~/Desktop/SRTR/Simulation/Github/kp-mods.Rdata")
source("~/Desktop/SRTR/Simulation/Github/generate_data.R") # generate_data function
source("~/Desktop/SRTR/getCI.R") # getCI function

sapply(a:b,function(i) {

# generate data
mydat.long = generate_data(i)

# first method:
## artificially censor when they deviate from the regime (no weights, not consistent for any parameter)
mydat.long = mydat.long %>% mutate(ACwt = ifelse(ld==1,0,1))
ACsurv = survfit(Surv(start,stop,death.timedep) ~ 1,data=mydat.long,
        type='kaplan-meier', weights=ACwt)
ACsum = summary(ACsurv,times=c(12,36,60)) # get estimated survival at the three time points
ACres = data.frame(time=c(12,36,60),
                   surv=ACsum[["surv"]],
                   lowCI=ACsum[["lower"]],
                   highCI=ACsum[["upper"]],
                   var.logsurv = (ACsum[["std.err"]] / ACsum[["surv"]])^2)
ACres$method = "AC"
# note: these are equivalent
# ACsum[["lower"]];ACsum[["upper"]]
# exp(log(ACres$surv) - qnorm(.975)*ACres$SE.logsurv);exp(log(ACres$surv) + qnorm(.975)*ACres$SE.logsurv)

# second method:
## artificially censor when they deviate from the regime AND
## weight by inverse probability of compliance
## (consistent for survival under representative intervention in study pop)

# subset to data up to tx
mydat.long.tx = mydat.long %>%
  filter(ld.lag==0,spk.lag==0)

# model probability of compliance
prob.ld.mod = glm(ld ~ ns(time,3) + dial_duration.timedep * pop,family=binomial,data=mydat.long.tx)
mydat.long[mydat.long$ld.lag==0 & mydat.long$spk.lag==0,"pred.prob.ld"] = prob.ld.mod$fitted.values
mydat.long[mydat.long$ld.lag!=0 | mydat.long$spk.lag!=0,"pred.prob.ld"] = 0
mydat.long = mutate(mydat.long,
                    ipACwt = ACwt * cumprod(1/(1-pred.prob.ld)))

# third method:
## artificially censor when they deviate from the regime AND
## weight by inverse probability of compliance AND
## adjust distribution of transplant rate to match target population

# assume ld is assigned first
# given that ld was not assigned, model prob of spk

# subset to data up to spk and before ld 
mydat.long.tx2 = mydat.long %>%
  filter(ld==0,spk.lag==0)

# model probability of spk
prob.spk.mod = glm(spk ~ ns(time,3) + dial_duration.timedep * pop,family=binomial,data=mydat.long.tx2)
mydat.long[mydat.long$ld==0 & mydat.long$spk.lag==0,"pred.prob.spk"] = prob.spk.mod$fitted.values
mydat.long[mydat.long$ld!=0 | mydat.long$spk.lag!=0,"pred.prob.spk"] = 0
mydat.long.pop1 = mydat.long
mydat.long.pop1$pop = 1
rate.probs.pop1 = predict(prob.spk.mod,subset(mydat.long.pop1,ld==0 & spk.lag==0),type="response")
mydat.long[mydat.long$ld==0 & mydat.long$spk.lag==0,"pred.prob.spk.pop1"] = rate.probs.pop1
mydat.long[mydat.long$ld!=0 | mydat.long$spk.lag!=0,"pred.prob.spk.pop1"] = 0
mydat.long = mutate(mydat.long,
                    prob.trt.received = ifelse(spk==1 & spk.lag==0,pred.prob.spk,1-pred.prob.spk),
                    prob.trt.received.pop1 = ifelse(spk==1 & spk.lag==0,pred.prob.spk.pop1,1-pred.prob.spk.pop1),
                    ip.rate.wt = ipACwt * cumprod(prob.trt.received.pop1/prob.trt.received))

# fourth method:
## artificially censor when they deviate from the regime AND
## weight by inverse probability of compliance AND
## adjust distribution of transplant rate to match target population AND
## adjust distribution of organ quality to match target population

# subset to spk transplant times
mydat.long.spk = mydat.long %>%
  filter(spk==1,spk.lag==0)

# quality depends only on pop; therefore, fit model for pop conditional on quality
kdri.mod = glm(pop ~ log.kdri,data=mydat.long.spk,family=binomial)

# marginal pop probs in this sample
marg.probs = prop.table(table(mydat.long[mydat.long$start==0,"pop"]))

# make quality weights
mydat.long.spk$prob.pop1.cond = kdri.mod$fitted.values
mydat.long.spk = mydat.long.spk %>% mutate(
  prob.pop.obs.cond = ifelse(pop==1,prob.pop1.cond,1-prob.pop1.cond),
  prob.pop.obs.marg = ifelse(pop==1,marg.probs["1"],marg.probs["0"]),
  qual.wt = (prob.pop1.cond / marg.probs["1"]) *
    (prob.pop.obs.marg / prob.pop.obs.cond)
)
mydat.long[mydat.long$spk==1 & mydat.long$spk.lag==0,"qual.wt"] = mydat.long.spk$qual.wt
mydat.long[mydat.long$spk==0 & mydat.long$spk.lag==0,"qual.wt"] = 1
mydat.long = mutate(mydat.long, qual.wt = na.locf(qual.wt,na.rm=F)) # fill in NAs
mydat.long = mutate(mydat.long,
                    ip.qual.wt = ip.rate.wt*qual.wt)

# fifth method:
## artificially censor when they deviate from the regime AND
## weight by inverse probability of compliance AND
## adjust distribution of transplant rate to match target population AND
## adjust distribution of KDRI to match target population AND
## adjust distribution of dialysis duration at baseline to match target population

## make covariate wts; condition on age and on_dial

# subset to listing time
mydat.bl = mydat.long %>% filter(time==0)

# fit model for pop conditional on age, on_dial
cov.mod = glm(pop ~ age + on_dial,data=mydat.bl,family=binomial)

# observed conditional probabilities
mydat.long[mydat.long$time==0,"cond.prob.pop1"] = cov.mod$fitted.values
mydat.long = mutate(mydat.long, cond.prob.pop1 = na.locf(cond.prob.pop1,na.rm=F))

# make covariate weights
mydat.long = mutate(mydat.long,
                    prob.pop.obs.cond = ifelse(pop==1,cond.prob.pop1,1-cond.prob.pop1),
                    prob.pop.obs.marg = ifelse(pop==1,marg.probs["1"],marg.probs["0"]),
                    cov.wt = (cond.prob.pop1 * prob.pop.obs.marg) /
                      (marg.probs["1"] * prob.pop.obs.cond),
                    ip.all.wt = ip.qual.wt * cov.wt)

## truncate all weights at 10
mydat.long = mutate(mydat.long,
                   ipACwt = ifelse(ipACwt > 10,10,ipACwt),
                   ip.rate.wt = ifelse(ip.rate.wt > 10,10,ip.rate.wt),
                   ip.qual.wt = ifelse(ip.qual.wt > 10,10,ip.qual.wt),
                   ip.all.wt = ifelse(ip.all.wt > 10,10,ip.all.wt))

# estimate weighted survival (representative intervention)
ipACsurv = survfit(Surv(start,stop,death.timedep) ~ 1,data=mydat.long,type='kaplan-meier', weights=ipACwt)
ipACsurv.CI = getCI(ipACsurv,mydat.long,"ipACwt")
ipACsurv.CI.sm = ipACsurv.CI[ipACsurv.CI$time%in%c(12,ipACsurv.CI$time[max(which(ipACsurv.CI$time <= 36))],max(ipACsurv.CI$time)),]
ipACsurv.CI.sm$method = "ipw-c"

# estimate weighted survival w/ rate weights only
ip.rate.surv = survfit(Surv(start,stop,death.timedep) ~ 1,data=mydat.long,type='kaplan-meier', weights=ip.rate.wt)
ip.rate.CI = getCI(ip.rate.surv,mydat.long,"ip.rate.wt")
ip.rate.CI.sm = ip.rate.CI[ip.rate.CI$time%in%c(12,ip.rate.CI$time[max(which(ip.rate.CI$time <= 36))],max(ip.rate.CI$time)),]
ip.rate.CI.sm$method = "ipw-rate"

# estimate weighted survival w/ rate and quality weights (should be consistent for study estimand)
ip.qual.surv = survfit(Surv(start,stop,death.timedep) ~ 1,data=mydat.long,type='kaplan-meier', weights=ip.qual.wt)
ip.qual.CI = getCI(ip.qual.surv,mydat.long,"ip.qual.wt")
ip.qual.CI.sm = ip.qual.CI[ip.qual.CI$time%in%c(12,ip.qual.CI$time[max(which(ip.qual.CI$time <= 36))],max(ip.qual.CI$time)),]
ip.qual.CI.sm$method = "ipw-qual"

# estimate weighted survival w/ all weights (should be consistent for target estimand)
ip.all.surv = survfit(Surv(start,stop,death.timedep) ~ 1,data=mydat.long,type='kaplan-meier', weights=ip.all.wt)
ip.all.CI = getCI(ip.all.surv,mydat.long,"ip.all.wt")
ip.all.CI.sm = ip.all.CI[ip.all.CI$time%in%c(12,ip.all.CI$time[max(which(ip.all.CI$time <= 36))],max(ip.all.CI$time)),]
ip.all.CI.sm$method="ipw-all"

# save results
res = rbind.data.frame(ACres,ipACsurv.CI.sm,ip.rate.CI.sm,ip.qual.CI.sm,ip.all.CI.sm)
res$iter = i

# write res to a text file
path = "~/Desktop/SRTR/Simulation/Github/"
write.table(res,file=paste0(path,"results.",i,".txt"),append=F,row.names=FALSE,col.names=T)
})
