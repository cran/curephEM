
basehaz.cureph = function(object)
{
  time=knots(object$baseline)
  hazard=object$baseline(time)

  return(data.frame(hazard=hazard,time=time))
}

survpred = function (object, newdata , time, center)
{
  UseMethod("survpred", object)
}

survpred.cureph = function(object, newdata , time, center = F)
{
  tt1 <- delete.response(object$terms.logistic)
  tt2 <- delete.response(object$terms.cox)

  if (missing(newdata) || is.null(newdata)) {
    A0 <- c(1,object$means$logistic * center)
    Z0 <- object$means$cox * center

    A <- model.matrix(tt1,object$data)
    Z <- model.matrix(tt2,object$data)[,-1]

    censored = object$Y[object$Y[,3]== -1,2]

    onesurv = T
    oneprob = F

  }else
  {
    censored = NULL

    A <- model.matrix(tt1,newdata)
    Z <- model.matrix(tt2,newdata)[,-1]

    oneprob = F
    if(nrow(A) == 1)
    {
      oneprob = T
      A0=A
    }else if(center){
      A0 = apply(A,2,mean,na.rm=T)
      binary.A = apply(A, 2,check.binary)
      A0[binary.A] = apply(A,2,min,na.rm=T)[binary.A]
      A0 = A0
    }else A0 = rep(0,ncol(A))
    A0[1]=1

    onesurv = F
    if(nrow(Z)==1)
    {
      onesurv = T
      Z0=Z
    }else if(center){
      Z0 = apply(Z,2,mean,na.rm=T)
      binary.Z = apply(Z, 2,check.binary)
      Z0[binary.Z] = apply(Z,2,min,na.rm=T)[binary.Z]
      Z0 = Z0
    }else Z0 = rep(0,ncol(Z))

  }

  if (missing(time) || is.null(time)) {
    time = knots(object$baseline)
  }else
    time = sort(unique(c(time,knots(object$baseline))))

  na.coef.logistic = is.na(object$coefficients$logistic)

  lor0 = A0[!na.coef.logistic] %*%
    object$coefficients$logistic[!na.coef.logistic]
  prob0 = exp(lor0)/(1+exp(lor0))
  if(!oneprob)
  {

  logistic.linear.predict =c(A[,!na.coef.logistic] %*%
                object$coefficients$logistic[!na.coef.logistic])
  names(logistic.linear.predict) = rownames(A)
  logistic.prob = exp(logistic.linear.predict)/(1+exp(logistic.linear.predict))
  }else{
    logistic.linear.predict = lor0
    logistic.prob = prob0
  }

  na.coef.cox = is.na(object$coefficients$cox)
  phaz0 = c(exp(Z0[!na.coef.cox]%*% object$coefficients$cox[!na.coef.cox]))
  cox.cumhaz = object$baseline(time)*phaz0

  cox.linear.predict=c(Z[,!na.coef.cox]%*% object$coefficients$cox[!na.coef.cox])
  if(!onesurv)
  {
    phazR = exp(cox.linear.predict) / phaz0
    surv.cox = exp(-outer(cox.cumhaz, phazR,'*'))
    surv.cureph = t(1-logistic.prob + logistic.prob * t(surv.cox))

    surv.cox = as.data.frame(cbind(time,surv.cox))
    surv.cureph = as.data.frame(cbind(time,surv.cureph))

    colnames(surv.cox) = colnames(surv.cureph) = c('Time', paste('Survival',rownames(A)))
  }else{
    surv.cox = exp(-cox.cumhaz)
    surv.cureph = 1-prob0 + prob0 * surv.cox

    surv.cox = as.data.frame(cbind(time,surv.cox))
    surv.cureph = as.data.frame(cbind(time,surv.cureph))

    colnames(surv.cox) = colnames(surv.cureph) = c('Time', "Survival")
  }

  surv.cox = rbind(surv.cox,c(object$end,rep(0,ncol(surv.cox)-1)))

out=list(logistic.linear.predict=logistic.linear.predict, logistic.prob=logistic.prob,
         cox.linear.predict=cox.linear.predict,
         cox.cumhaz=cox.cumhaz,surv.cox=surv.cox,surv.cureph=surv.cureph,
         origin = object$origin, end = object$end,censored = censored)
class(out)='survpred.cureph'
  return(out)
}

plot.survpred.cureph = function(x, pooled = T, censor = x$censored,...)
{
  if(pooled)
    surv = x$surv.cureph
  else
    surv = x$surv.cox

  plot(0,col=0,ylim = c(0,1),xlim = c(x$origin,x$end),
       ylab='Survival',xlab = 'Time')

  abline(0,0,col='grey')
  abline(1,0,col='grey')
  abline(v=x$end,lty=2)
  for(i in 2:ncol(surv))
  {
    tmpstep = stepfun(surv[,1],c(1,surv[,i]))
    lines(tmpstep,verticals= T, do.points = F,xlim = c(x$origin,x$end))
  }

  if(any(as.logical(censor)))
  {
    points(censor,tmpstep(censor),pch=4)
  }
}
