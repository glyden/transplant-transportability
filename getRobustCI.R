## compute a robust variance for the log of time-varying weighted Kaplan-Meier survival
# note: dat must have an id variable called "PERS_ID" and an outcome variable called "death.timedep," which cannot have NAs
getVar = function(fit,dat,wt) {
  # compute hazards
  sum = summary(fit)
  haz = sum$n.event / sum$n.risk
  n = sum(dat$time==0)
  n.times = length(sum$n.risk)
  
  ## make n*A matrix
  # make diagonal matrix of weighted sum of number at risk
  A = diag(sum$n.risk)
  # add extra columns for log survival (1 per time)
  A = cbind(A,matrix(0,nrow=n.times,ncol=n.times))
  # add bottom rows
  inv.haz = 1/(1-haz)
  survA = sapply(1:n.times,function(t) {
    c = rep(0,n.times)
    c[t] = 1
    vect = c(inv.haz[1:t],rep(0,n.times-t),c)
    return(vect)
  })
  A = rbind(A,n * t(survA))
  
  ## make n*B matrix
  # get psi for hazards
  ids = unique(dat$PERS_ID)
  psi_main = sapply(sum$time,function(x) {
    psi = rep(0,length(ids))
    names(psi) = ids
    psi[match(dat$PERS_ID[dat$stop==x],names(psi))] = 
      dat[[wt]][dat$stop==x] * (dat$death.timedep[dat$stop==x] - haz[which(sum$time==x)])
    return(psi)
  })
  # add the estimating equation for log survival to the stack (evaluated at estimates therefore zero)
  psi_all = cbind(psi_main,matrix(0,nrow=dim(psi_main)[1],ncol=n.times))
  # for each person, compute psi %*% t(psi), then sum across all people
  mats = lapply(1:nrow(psi_all),function(row) psi_all[row,] %*% t(psi_all[row,])) # should give us a list of matrices
  B = Reduce('+', mats)
  
  ## compute estimated variance of log survival at survtime
  varcov = solve(A) %*% B %*% t(solve(A))
  var.logsurv = (diag(varcov))[(n.times+1):(n.times*2)]
  return(var.logsurv)
}

## get 95% CIs based on the log transform
getCI = function(fit,dat,wt,max.time=60) {
  # get variance of log survival
  var.logsurv = getVar(fit,dat,wt)
  
  # compute CI at each time point
  sum = summary(fit)
  times = c(0,sum$time); survs=c(1,sum$surv)
  lowCI = c(1,exp(log(sum$surv) - qnorm(.975)*sqrt(var.logsurv)))
  highCI = c(1,pmin(exp(log(sum$surv) + qnorm(.975)*sqrt(var.logsurv)),1))
  
  # return output
  output = data.frame(time=times,surv=survs,lowCI=lowCI,highCI=highCI,var.logsurv=c(0,var.logsurv))
  output = mutate(output,lowCI = ifelse(surv==0,NA,lowCI),
                  highCI = ifelse(surv==0,NA,highCI))
  if (max(output$time) < max.time) {
    last.row = output[dim(output)[1],]
    new.rows = as.data.frame(lapply(last.row, rep, max.time - max(output$time)))
    new.rows$time = seq(max(output$time)+1,max.time)
    output = rbind.data.frame(output,new.rows)
  }
  return(output)
}
