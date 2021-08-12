### create getCI function ###

# function to get CIs for K-M survival estimator with time-varying weights
getCI = function(fit,dat,wt) {
  # compute variance using my formula (treating weights as fixed)
  sum = summary(fit)
  m.inv1 = sapply(sum$time,function(x) sum((dat[[wt]][dat$stop==x])^2))
  m.inv = m.inv1 / (sum$n.risk^2)
  haz = sum$n.event / sum$n.risk
  var.logsurv = cumsum(m.inv * (haz/(1-haz)))
  
  # compute CI at each time point
  times = c(0,sum$time); survs=c(1,sum$surv)
  lowCI = c(1,exp(log(sum$surv) - 1.96*sqrt(var.logsurv)))
  highCI = c(1,pmin(exp(log(sum$surv) + 1.96*sqrt(var.logsurv)),1))
  
  # return output
  output = data.frame(time=times,surv=survs,lowCI=lowCI,highCI=highCI,var.logsurv=c(0,var.logsurv))
  output = mutate(output,lowCI = ifelse(surv==0,NA,lowCI),
                  highCI = ifelse(surv==0,NA,highCI))
  return(output)
}
