# generate data for a simulation study with two populations with equal follow-up (five years)
# but different distributions of time on dialysis (specifically, proportion on dialysis at listing), different ages,
# different wait times to deceased-donor and living-donor transplant, and same survival conditional on dialysis time, 
# wait time to transplant, and organ quality

# in population 1:
# older ages; more on dialysis at listing
# faster rate of deceased-donor transplant
# slower rate of living-donor transplant
# better quality (lower KDRI)

generate_data = function(seed) {
  set.seed(seed)
  
  # set sample size
  n = 2000
  
  # assign age, whether people are on dialysis
  mydat = data.frame(age = sample(c(age_dist_worse,age_dist_better),n,replace=T),
                     on_dial = rbinom(n,1,.65))
  
  # get population probability conditional on age and whether people are on dialysis
  pop.x = as.matrix(cbind.data.frame(Intercept=rep(1,nrow(mydat)),mydat[,c("age","on_dial")]))
  expit = function(x) {
    return(exp(x) / (1+exp(x)))
  }
  mydat$prob.pop1 = expit(pop.x %*% pop_mod.coef)
  
  # randomly assign to a population
  mydat$pop = rbinom(nrow(mydat),size=1,prob=mydat$prob.pop1)
  num0 = unname(table(mydat$pop)["0"])
  
  ## generate months on dialysis
  
  # among those on dialysis, assign how long they've been on (in days)
  mydat$dial_duration = 0
  mydat$dial_duration[mydat$on_dial==1] = sample(dial_time,sum(mydat$on_dial==1),replace=T)
  
  # among those not on dialysis, assign how long they wait to start
  mydat$dial_start = NA
  mydat$dial_start[mydat$on_dial==0] = sample(time_to_dial_start,sum(mydat$on_dial==0),replace=T)
  
  # make long data
  mydat$id = 1:nrow(mydat)
  mydat.long = expandRows(mydat,count=60,count.is.col=F,drop=F)
  mydat.long$time = rep(0:59,n)
  ns.time = ns(mydat.long$time,knots=c(5,17),Boundary.knots=c(0,212))
  colnames(ns.time) = paste0("ns.time",1:3)
  mydat.long = cbind.data.frame(mydat.long,ns.time) # get ns(time,3) for wl_surv_mod, spk_mod, ld_mod
  mydat.long$start = mydat.long$time
  mydat.long$stop = mydat.long$start+1
  mydat.long = mutate(mydat.long,
                      dial_duration.mths = ifelse(dial_duration==0 & !is.na(dial_start),
                                                  -dial_start/30,dial_duration/30),
                      dial_duration.timedep = ifelse(dial_duration.mths==0,0,
                                                     pmax(0,time+dial_duration.mths)))
  
  # get discrete time hazards at each time conditional on dialysis duration
  wl_surv.x = as.matrix(cbind.data.frame(Intercept=rep(1,nrow(mydat.long)),
                                         mydat.long[,c("ns.time1","ns.time2","ns.time3","age",
                                                       "dial_duration.timedep")]))
  mydat.long$prob.wl.dth = expit(wl_surv.x %*% wl_surv_mod.coef)
  
  # assign waitlist failure
  mydat.long$wl.dth.assign = rbinom(nrow(mydat.long),size=1,prob=mydat.long$prob.wl.dth)
  mydat.long = mydat.long %>% group_by(id) %>% mutate(wl.dth = pmin(1,cumsum(wl.dth.assign)))
  mydat.long = mutate(mydat.long, wl.dth.lag = ifelse(row_number()==1,0,lag(wl.dth)))
  
  ## assign time to dd and ld transplant
  
  # assign ld transplant
  ld.x = as.matrix(cbind.data.frame(Intercept=rep(1,nrow(mydat.long)),
                                    mydat.long[,c("ns.time1","ns.time2","ns.time3",
                                                  "dial_duration.timedep","pop")],
                                    int1=mydat.long$dial_duration.timedep*mydat.long$pop))
  mydat.long$prob.ld = expit(ld.x %*% ld_mod.coef)
  mydat.long$prob.ld = pmin(.25,4*mydat.long$prob.ld)
  mydat.long$ld.assign = rbinom(nrow(mydat.long),size=1,prob=mydat.long$prob.ld)
  mydat.long = mydat.long %>%
    mutate(ld = pmin(1,cumsum(ld.assign)))
  
  # assign dd transplant
  spk.x = as.matrix(cbind.data.frame(Intercept=rep(1,nrow(mydat.long)),
                                     mydat.long[,c("ns.time1","ns.time2","ns.time3",
                                                   "dial_duration.timedep","pop")],
                                     int=mydat.long$dial_duration.timedep*mydat.long$pop))
  mydat.long$prob.spk = expit(spk.x %*% spk_mod.coef)
  mydat.long$spk.assign = rbinom(nrow(mydat.long),size=1,prob=mydat.long$prob.spk)
  mydat.long = mydat.long %>% 
    mutate(spk = pmin(1,cumsum(spk.assign)))
  
  # remove rows after ld transplant, if not previously transplanted with dd (don't use these)
  mydat.long = mydat.long %>% mutate(ld = ifelse(ld==1 & spk==1 & cumsum(spk)>cumsum(ld),0,ld),
                                     ld.lag = ifelse(row_number()==1,0,lag(ld))) %>%
    filter(ld.lag==0) %>%
    mutate(spk = ifelse(ld==1,0,spk)) # assume ld is assigned first, then spk
  
  # assign organ quality to dd transplants
  mydat.long$log.kdri = NA
  mydat.long = mutate(mydat.long,spk.lag = ifelse(row_number()==1,0,lag(spk)))
  mydat.long$log.kdri[mydat.long$spk==1 & mydat.long$spk.lag==0 & mydat.long$pop==1] = 
    rnorm(n=sum(mydat.long$spk==1 & mydat.long$spk.lag==0 & mydat.long$pop==1),
          mean=-.1721,sd=.1282)
  mydat.long$log.kdri[mydat.long$spk==1 & mydat.long$spk.lag==0 & mydat.long$pop==0] = 
    rnorm(n=sum(mydat.long$spk==1 & mydat.long$spk.lag==0 & mydat.long$pop==0),
          mean=-0.0013,sd=.1282)
  mydat.long = mutate(mydat.long, log.kdri = na.locf(log.kdri,na.rm=F)) # fill in NAs

  # exponentiate log KDRI
  mydat.long = mutate(mydat.long,kdri = exp(log.kdri))

  ## assign post-dd-transplant survival by drawing from covariate-specific survival curves
  # (inverse transform sampling)
  mydat.long.spk = mydat.long %>% filter(spk==1,spk.lag==0,wl.dth.lag==0) %>% rename(spk_KDRI=kdri)
  mydat.long.spk$posttx_unif = runif(nrow(mydat.long.spk))
  mydat.long.spk$age = mydat.long.spk$age + (mydat.long.spk$start/12)
  posttx.survs = lapply(1:nrow(mydat.long.spk),function(i) {
    unif.prob = mydat.long.spk$posttx_unif[i]
    x = as.numeric(mydat.long.spk[i,c("age","dial_duration.timedep","spk_KDRI")])
    risk.score = exp(x %*% posttx.coef)
    surv_curve = exp(-(as.vector(risk.score) * posttx.bh$hazard))
    # surv_curve2 = survfit(posttx.mod,newdata=mydat.long.spk[i,]) # same
    cdf = 1-surv_curve
    if (all(cdf<=unif.prob)) {
      posttx_surv = max(posttx.bh$time)
      posttx_cens = 1
    } else {
      posttx_surv = posttx.bh$time[which(cdf > unif.prob)[1]]
      posttx_cens = 0
    }
    return(data.frame(posttx_surv=posttx_surv,posttx_cens=posttx_cens))
  })
  posttx.survs = do.call(rbind.data.frame,posttx.survs)
  posttx.survs.id = cbind.data.frame(id=mydat.long.spk$id,posttx.survs)
  mydat.long = merge(mydat.long,posttx.survs.id,by="id",all.x=T) %>% arrange(id,time) # no longer grouped
  mydat.long = mydat.long %>% group_by(id) %>%
    mutate(posttx_dth = ifelse(spk==1 & cumsum(spk.lag)>=posttx_surv & posttx_cens==0,1,0))
  mydat.long = mutate(mydat.long,
                      death.timedep = ifelse(ld==1,NA,
                                             ifelse(cumsum(wl.dth) > cumsum(spk),1,
                                                    ifelse(posttx_dth==1,1,0))))
  
  # remove rows after death
  mydat.long = mutate(mydat.long,death.timedep.lag = ifelse(row_number()==1,0,lag(death.timedep))) %>%
    filter(death.timedep.lag==0)
  
  # remove variables
  mydat.long = mydat.long %>% select(-dial_start,-contains("ns"),-dial_duration.mths,
                                     -prob.wl.dth,-wl.dth.assign,-prob.ld,-ld.assign,-prob.spk,-spk.assign,
                                     -posttx_surv,-posttx_cens)
  return(mydat.long)
}