
cureph.simgen = function()
{
  dat=CureCox.simgen()
  sim.cureph.data = data.frame(row.names = 1:200)
  sim.cureph.data$time=dat$survn[,1]
  sim.cureph.data$time2=dat$survn[,2]
  sim.cureph.data$event=dat$survn[,3]
  sim.cureph.data$Z1 = dat$An[,1]
  sim.cureph.data$Z2 = factor(dat$An[,2])
  Z3=factor(c('A','B','C')[sample(3,200,replace = T)])
  sim.cureph.data$Z3=Z3
  sim.cureph.data$Z4 = factor((Z3 == 'A') + (dat$An[,2]==1))

  attr(sim.cureph.data,'true.coef')=list(
    a = c(1,-0.63,1,rep(0,4)),b=c(-0.2,0.3,rep(0,4))
  )
  attr(sim.cureph.data,'true.baseline.surv') =
    function(x) pmax(pmin(1-x/20,1),0)

  return(sim.cureph.data)
}

inv.T.=function(c,tao)
{
  u=runif(length(c))
  t=-log(u)/c
  return(inv.L0(t,tao))
}

inv.L0=function(t,tao)
{
  return(tao*(1-exp(-t)))
}

cov.gen.=function(N,pA=2,pZ=2)
{
  Z=cbind(rnorm(N,4),rbinom(N,1,0.3))

  return(list(A=Z,Z=Z))
}

Q2C2.gen=function(N,tao)
{
  Q=runif(N,max=13)#*rbinom(N,1,0.95)
  C=runif(N,min=16,max=22)
  return(list(Q=Q,C=C))
}

CureCox.simgen=function(n=200,pA,pZ,tao=20,N,a=c(1,-0.63,1),b=c(-0.2,0.3),
                        inv.T=inv.T.,QC.gen=Q2C2.gen,
                        cov.gen=cov.gen.)
{
  if(missing(pA))
    pA=length(a)-1

  if(missing(pZ))
    pZ=length(b)

  if(is.null(a))
    a=rnorm(pA)

  if(is.null(b))
    b=rnorm(pZ)

  if(missing(N))
    N=2*n

  survn=yn=An=Zn=NULL

  Ncount=Qcount=Ccount=0

  while(n>0)
  {
    tmp.cov=cov.gen(N,pA,pZ)
    A=tmp.cov$A # generate glm covariates
    Z=tmp.cov$Z # generate coxph covariates

    or=exp(a[1]+A%*%a[-1]) # glm odds ratio
    pn=or/(1+or)      # uncure probability

    ph=exp(Z%*%b)     # cox hazard ratio

    Tn=rep(Inf,N) # initialize event time

    Yn=rbinom(N,1,pn) # generate cure portion
    Ncount=Ncount+sum(Yn)
    # Q and C might cause trouble with S(Q) or S(C) being zero!
    QC=QC.gen(N,tao)
    Qn=QC$Q # uniform truncation
    #Qn=rep(0,N)
    Cn=QC$C  # independent residual censoring
    #Cn=rep(tao+1,N)
    Tn[Yn==1]=inv.T(ph[Yn==1],tao) # uncured event times

    event= Tn<=Cn   # observe censoring
    Ccount=Ccount+sum(event)
    Y=-1+2*event+(Yn==0)*(Cn>=tao) # observe Y

    # survival data
    time1=Qn
    time2=pmin(Tn,Cn)

    # truncation
    obs= Tn>Qn
    Qcount=Qcount+sum(!obs)
    m=min(sum(obs),n)
    n=n-m
    obs= which(obs)[1:m]

    # attach
    survn=Surv(c(survn[,1],time1[obs]),
               c(survn[,2],time2[obs]),
               c(survn[,3],event[obs]))
    yn=c(yn,Y[obs])
    An=rbind(An,matrix(A[obs,],m,pA))
    Zn=rbind(Zn,matrix(Z[obs,],m,pZ))
  }

  Qrate=Qcount/Ncount
  Crate=1-Ccount/Ncount

  return(list(survn=survn,yn=yn,An=An,Zn=Zn,
              Qrate=Qrate,Crate=Crate))
}
