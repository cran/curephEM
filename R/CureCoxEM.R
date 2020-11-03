
# Main Function
#################################################################
# survn: Survival object of "interval" type
# yn: observed cure vector (1-uncured, 0-cured, -1-censored)
# X: covariates for glm
# Z: covariates for coxph
# a0,b0,coxfit0: initial values


CureCox.EM=function(Y, A, Z, init, control, var.method)
{
  a=init[[1]]
  b=init[[2]]
  nvar=1+ncol(A)+ncol(Z)
  nobs=nrow(A)
  flag1=rep(T,ncol(A))
  if(length(flag1)>0)
  {
    tmpA=cbind(1,A)
    cov.chol=suppressWarnings(chol(t(tmpA)%*% tmpA,pivot = T, tol = control$toler.chol))

    while(attr(cov.chol,'rank') < (1+sum(flag1)))
    {
      for(j in sort(which(flag1), decreasing = T))
      {
        tmp.keep=flag1
        tmp.keep[j]=F
        tmpA=cbind(1,A[,tmp.keep])

        tmp.chol=suppressWarnings(chol(t(tmpA)%*% tmpA,pivot = T, tol = control$toler.chol))

        if(attr(tmp.chol,'rank') >= attr(cov.chol,'rank'))
        {
          flag1=tmp.keep
          cov.chol=tmp.chol

          break
        }
      }
    }
  }

  flag2=rep(T,ncol(Z))
  if(length(flag2)>0)
  {
    tmpZ=cbind(1,Z)
    cov.chol=suppressWarnings(chol(t(tmpZ)%*% tmpZ,pivot = T, tol = control$toler.chol))

    while(attr(cov.chol,'rank') < (1+sum(flag2)))
    {
      for(j in sort(which(flag2), decreasing = T))
      {
        tmp.keep=flag2
        tmp.keep[j]=F
        tmpZ=cbind(1,Z[,tmp.keep])

        tmp.chol=suppressWarnings(chol(t(tmpZ)%*% tmpZ,pivot = T, tol = control$toler.chol))

        if(attr(tmp.chol,'rank') >= attr(cov.chol,'rank'))
        {
          flag2=tmp.keep
          cov.chol=tmp.chol

          break
        }
      }
    }
  }

  yn=Y[,3]
  survn=survival::Surv(Y[,1],Y[,2],as.numeric(yn==1))
  tmp = CureCox.EM.fit(survn, yn,A[,flag1],Z[,flag2],control$iter.max,control$eps,
              a[c(TRUE,flag1)],b[flag2],init[[3]])
  a[c(TRUE,flag1)] = tmp$a
  a[c(FALSE,!flag1)] = NA
  b[!flag2]=NA
  b[flag2]=tmp$b
  baseline=tmp$Lam0

  coef=c(a,b)
  flag=c(TRUE,flag1,flag2)
  var=matrix(0,nvar,nvar)

  var[flag,flag]=var.method(survn,yn,matrix(A[,flag1],nobs),
                            matrix(Z[,flag2],nobs),
                            a[c(TRUE,flag1)],b[flag2],baseline)

  iter=tmp$iter
  loglik=CureCox.EM.loglik(tmp,Y,Y[,3],matrix(A[,flag1],nobs),
                           matrix(Z[,flag2],nobs))

  #    concordance <- survConcordance.fit(y, lp, strata, weights)
  list(coefficients = coef, var = var, loglik = loglik,baseline = tmp$Lam0,
       #         score = agfit$sctest,
       iter = iter,
       #          linear.predictors = as.vector(lp),
       #         residuals = resid, means = agfit$means, concordance = concordance,
       method = "curph.EM")

}

CureCox.EM.fit=function(survn,yn,A=NULL,Z=NULL,nstep=200,thres=1e-8,
                 a0=NULL,b0=NULL,Lam0=NULL)
{
  N=length(yn)  # sample size

  if(is.null(A))
    A=matrix(0,N,0)
  else
    A=matrix(A,N)

  pA=ncol(A)
  if(is.null(Z))  # two models share same covariates
    Z=A
  Z=matrix(Z,N)
  pZ=ncol(Z)

  if(is.null(a0))   # initialize a0 (with intercept)
    a0=coef(glm(yn~A,subset= yn!=-1)) # start all 0 by default

  if(pZ==0)
    coxfit0=survival::coxph(survn~1,subset= yn==1,ties="breslow") # fit a naive coxph
  else
    coxfit0=survival::coxph(survn~Z,subset= yn==1,ties="breslow") # fit a naive coxph

  if(is.null(b0))   # initialize b0
    b0=coef(coxfit0)  # start all 0 by default

  tk=survival::survfit(coxfit0)$time # unique event times
  if(is.null(Lam0)) # initialize baseline hazard
  {
    tmpBaseH=survival::basehaz(coxfit0,centered=F)
    Lam0=stepfun( tk,c(0,tmpBaseH$hazard[is.element(tmpBaseH$time,tk)]))
  }
#  Lam0=-log(1-tk/(max(tk)+1))
  K=length(tk)        # number of unique event times
  dn= yn==-1   # censoring indicator

  for(i in 1:nstep)
  {
    #print(i)

    or=exp(a0[1]+A%*%a0[-1]) # odds ratio
    pn=or/(1+or)           # probability
    if(pZ==0)
      ph=rep(1,N)
    else
      ph=exp(Z%*%b0)    # hazard ratio

    # E-step
    tempE=CureCox.E(survn,yn,pn,ph,Lam0,tk)

    # M-step
    tempM=CureCox.M(tk,survn[dn,2],dn,A,Z,tempE$wy1,tempE$wy0,tempE$wfnk,tempE$wSn)

    # check convergence
    oldhat=c(a0,b0,Lam0(tk))

    # update parameters
    a0=tempM$a0
    b0=tempM$b0
    curefit0=tempM$curefit0
    coxfit0=tempM$coxfit0

    wSi=rep(0,N)
    wSi[dn]=tempE$wSn
    if(pZ==0)
      ph=rep(1,N)
    else
      ph=exp(Z%*%b0)

    Lam0=CureCox.basehaz(survn[,2],tk,tempE$wfnk,tempE$wfk,wSi,ph)

      newhat=c(a0,b0,Lam0(tk))
      move=max(abs((oldhat-newhat)))
#      print(c(i,move,a0,b0))
      if(move<thres)
      {
        message(paste("Converge at step ",i))
        break
      }
  }

  if(move>thres)
    warning(paste("Fail to converge after step ",nstep))

  lam0=diff(Lam0(c(0,tk)))

  return(list(a=coef(curefit0),b=coef(coxfit0),Lam0=Lam0,lam0=lam0))
}

# E-step returns weights
# survn: observed survival data
# yn: observed cure vector (1-uncured, 0-cured, -1-censored)
# pn: estimated cure rate vector
# fnk: estimated mass at event times N*K matrix
# tk: observed event time
CureCox.E=function(survn,yn,pn,ph,Lam0,tk)
{
  N=nrow(survn) # Sample size
  dn= yn==-1   # censor indicator
  cN=sum(dn)   # censor sample size
  K=length(tk) # mass points

  # I(qn>tk)
  I.tk.qn= (rep(1,K)%*%t(survn[,1])) > (tk%*%t(rep(1,N)))


  Lnk=Lam0(tk)%*%t(ph)
  Snk=exp(-Lnk)
  lam0=diff(Lam0(c(0,tk)))
  fnk=outer(lam0,c(ph))*Snk   # event density matrix

  Fqn=apply(fnk*I.tk.qn,2,sum)   # P(Ti<qi)
  mn=pn*Fqn/(1-pn*Fqn)                # ghost sample size
#  hist(mn)

  if(cN>0)
  {
    Sncn=exp(-Lam0(survn[dn,2])*ph[dn])  # Sn(cn)
  }
  Fnqn.inv=rep(1,K)%*%t(1/Fqn)   # 1/Fn(qn)
  Fnqn.inv[is.infinite(Fnqn.inv)]=0 # remove Inf

  if(cN>0)
    Eyn=pn[dn]*Sncn/(pn[dn]*Sncn+1-pn[dn])

  # Weight for log p
  wy1=(yn==1)  # observe y=1
  if(cN>0) wy1[dn]=Eyn # censored y=-1
  wy1=wy1+mn # truncated ghost samples

  # Weight for log(1-p)
  wy0=(yn==0) # observe y=0
  if(cN>0) wy0[dn]=1-Eyn # censored y=-1

  # Weight for log fn(tk)
  wfnk=matrix(0,K,N) # initialization
  wfnk[,!dn]=       # observed event times
    (tk %*% t(rep(1,N-cN)))==( rep(1,K)%*% t(survn[!dn,2]))

  wfnk=wfnk+I.tk.qn*fnk*Fnqn.inv*(rep(1,K)%*%t(mn))  # truncated ghost samples
  wfk=apply(wfnk,1,sum)
  if(cN>0) wSn=Eyn else wSn=NULL  # weight for Sn(Cn)

  return(list(wy1=wy1,wy0=wy0,wfnk=wfnk,wSn=wSn,wfk=wfk))

}

# M-step update parameters
# tk: observed event times
# A: glm covariates
# Z: coxph covariates
# wy1,wy0: glm weights
# wfnk: coxph weights
CureCox.M=function(tk,tc,dn,A,Z,wy1,wy0,wfnk,wSn)
{
  N=length(wy1) # sample size
  K=length(tk)  # mass points

  # select non-zero weights
  idy1= wy1!=0  #
  idy0= wy0!=0  #

  # create data for glm
  wy=c(wy1[idy1],wy0[idy0])    # wy11,wy12,...,wy1n,wy01,...,wy0n
  y=rep(1:0,c(sum(idy1),sum(idy0))) # 1,1,...,1,1,0,0,...,0,0
  A=matrix(A[c(which(idy1),which(idy0)),],ncol=ncol(A))     # x1,x2,...,xn,x1,x2,...,xn

  # fit glm
  if(ncol(A)==0)
    curefit=suppressWarnings(glm(y~1,family=binomial,weights=wy))
  else
    curefit=suppressWarnings(glm(y~A,family=binomial,weights=wy))
  a0=coef(curefit)
  or=predict(curefit)

  # select non-zero weights
  idcoxf= wfnk != 0
  idcoxS= wSn != 0

  # create data for coxph
  en=sum(idcoxf)
  cn=sum(idcoxS)

  s=survival::Surv(
      time=c(tk[unlist(apply(idcoxf,2,which))], # t1,t2,...,tK,t1,...,tK,t1,...
        tc[which(idcoxS)]), event=rep(c(T,F),c(en,cn)))    # censoring time
  Z=matrix(Z[c(rep(1:N,apply(idcoxf,2,sum)),  # Z1,Z1,...,Z1,Z2,...,Z2,Z3,...
        which(dn)[idcoxS]),],ncol=ncol(Z))     # Z for censored obs
  wcox=c(wfnk[idcoxf],         # wfnk11,wfnk12,...,wfnk1K,wfnk21,.....
         wSn[idcoxS])     # weight 1 for censored obs

  # fit coxph
  if(ncol(Z)==0)
    coxfit0=survival::coxph(s~1,weights=wcox)
  else
    coxfit0=survival::coxph(s~Z,weights=wcox)

  b0=coef(coxfit0)  # get coefficients


  return(list(a0=a0,b0=b0,coxfit0=coxfit0,curefit0=curefit))
}

CureCox.basehaz=function(X,tk,wfnk,wfk,wSi,ph)
{
  N=length(X)
  K=length(tk)

#  wfk=apply(wfnk,1,sum)
  wfij= rep(1,K)%*%t(apply(wfnk,2,sum))-apply(rbind(0,wfnk[-K,]),2,cumsum)

  tmp=rep(1,K) %*% t(wSi)
  I.Xi.tk= (rep(1,K)%*%t(X)) >= (tk%*%t(rep(1,N)))
  tmp=tmp * I.Xi.tk+wfij

  tmp=tmp %*% diag(c(ph))

  lhat=wfk/apply(tmp,1,sum)

  return(stepfun(tk,c(0,cumsum(lhat))))
}

CureCox.EM.loglik=function(fit,survn,yn,Xn=NULL,Zn=NULL,a=fit$a,b=fit$b,Lam0=fit$Lam0)
{
  tk=knots(fit$Lam0)
  dt=min(diff(tk))/2
  Xn=cbind(1,Xn)

  Xa=Xn %*% a
  if(is.null(Zn))
    Zb=rep(0,length(yn))
  else
    Zb=Zn %*% b
  LamX=Lam0(survn[,2]) * exp(Zb)
  LamQ=Lam0(survn[,1]) * exp(Zb)
  logdLamX=ifelse(yn==1,log(Lam0(survn[,2])-Lam0(survn[,2]-dt)),0)

  out= (yn==1)*(Xa+Zb + logdLamX
                - LamX)+
    (yn==-1)*log(1+exp(Xa-LamX))-
    log(1+exp(Xa-LamQ))

  return(sum(out))
}
