CureCox.Louis.varMatrix=function(survn,yn,A,Z,a,b,Lam0)
{
  gradn=CureCox.Louis.Egradn(survn,yn,A,Z,a,b,Lam0)

  Hess=CureCox.Louis.Ehess(survn,yn,A,Z,a,b,Lam0)

  grad.sq=CureCox.Louis.Egradsq(survn,yn,A,Z,a,b,Lam0)

  #  print(max(abs(apply(gradn,2,sum))))
  #  print(min(eigen(-Hess[1:5,1:5])$value))
  #  print(min(eigen(-Hess)$value))

  N=nrow(survn) # Sample size
  ker=rep(1,N)%*%t(rep(1,N))
  diag(ker)=0
  Louis=-Hess-t(gradn)%*%ker%*%gradn-grad.sq


  #out=t(gradn)%*%ker%*%gradn+grad.sq

  out=as.matrix(solve(Louis)[1:(length(a)+length(b)),1:(length(a)+length(b))])
  return(out)
}

CureCox.Louis.var=function(survn,yn,A,Z,a,b,Lam0)
{
  if(any(is.na(a))|any(is.na(b))|is.null(Lam0))
    return(NA*c(a,b))
  else
    return(diag(CureCox.Louis.varMatrix(survn,yn,A,Z,a,b,Lam0)))
}

# Hessian of log-likelihood
CureCox.Louis.Ehess=function(survn,yn,A,Z,a,b,Lam0)
{
  tk=knots(Lam0)
  K=length(tk)        # number of unique event times
  nk=table(survn[yn==1,2]) # counting ties, usually all 1
  dn= yn==-1   # censoring indicator
  N=nrow(survn) # Sample size
  cN=sum(dn)   # censor sample size

  A=cbind(1,A)

  # Logisitc Probability
  eor=c(exp(A%*%a))
  pn=eor/(eor+1)

  # Hazard Proportion
  phaz=c(exp(Z%*%b))

  # Hazard and Density
  Lnk=phaz%*%t(Lam0(tk))
  Snk=exp(-Lnk)
  lamk=diff(c(0,Lam0(tk)))
  col.lamk=(rep(1,N)%*%t(lamk))

  # Indicators
  I.Qi.tk= (survn[,1]%*%t(rep(1,K))) > (rep(1,N)%*%t(tk))
  I.Xi.tk= (survn[,2]%*%t(rep(1,K))) >= (rep(1,N)%*%t(tk))
  I.Xieqtk= (survn[,2]%*%t(rep(1,K))) == (rep(1,N)%*%t(tk))

  # Expectation on latent variables
  Eyn=yn
  tmp=eor[dn]*exp(-phaz[dn]*Lam0(survn[dn,2]))
  Eyn[dn]=tmp/(1+tmp)

  fnk=outer(c(phaz),lamk)*Snk
  Fqn=apply(fnk*I.Qi.tk,1,sum)
  tmp=eor*Fqn
  Emn=tmp/(1+eor-tmp)

  # Ghost Probability
  PTitk=fnk*I.Qi.tk
  PTitk=PTitk/Fqn
  PTitk[is.nan(PTitk)]=0
  PTijtk=t(1-apply( cbind(0,PTitk[,-K]) , 1, cumsum))

  # Hessian for aa
  waa=c((1+Emn)*pn*(1-pn))
  Haa=-t(A)%*% diag(waa) %*% A

  # Hessian for bb
  wbb=+Eyn*phaz*Lam0(survn[,2])+Emn*apply(PTitk*Lnk,1,sum)
  Hbb=-t(Z) %*% diag(wbb)%*% Z

  # Hessian for blamb

  wblam=(I.Xi.tk*Eyn+PTijtk*Emn)*phaz
  Hblam=-t(Z)%*%wblam

  # Hessian for lamlam
  Hll=diag(-(nk+apply(PTitk*Emn,2,sum))/lamk^2)

  temp=rbind(cbind(Hbb,Hblam),cbind(t(Hblam),Hll))

  return(Matrix::bdiag(Haa,temp))
}

# Gradient of log-likelihood
CureCox.Louis.Egradn=function(survn,yn,A,Z,a,b,Lam0)
{
  tk=knots(Lam0)
  K=length(tk)        # number of unique event times
  nk=table(survn[yn==1,2]) # counting ties, usually all 1
  dn= yn==-1   # censoring indicator
  N=nrow(survn) # Sample size
  cN=sum(dn)   # censor sample size

  A=cbind(1,A)

  # Logisitc Probability
  eor=c(exp(A%*%a))
  pn=eor/(eor+1)

  # Hazard Proportion
  phaz=c(exp(Z%*%b))

  # Hazard and Density
  Lnk=phaz%*%t(Lam0(tk))
  Snk=exp(-Lnk)
  lamk=diff(c(0,Lam0(tk)))
  col.lamk=(rep(1,N)%*%t(lamk))

  # Indicators
  I.Qi.tk= (survn[,1]%*%t(rep(1,K))) > (rep(1,N)%*%t(tk))
  I.Xi.tk= (survn[,2]%*%t(rep(1,K))) >= (rep(1,N)%*%t(tk))
  I.Xieqtk= (survn[,2]%*%t(rep(1,K))) == (rep(1,N)%*%t(tk))

  # Expectation on latent variables
  Eyn=yn
  tmp=eor[dn]*exp(-phaz[dn]*Lam0(survn[dn,2]))
  Eyn[dn]=tmp/(1+tmp)

  fnk=outer(c(phaz),lamk)*Snk
  Fqn=apply(fnk*I.Qi.tk,1,sum)
  tmp=eor*Fqn
  Emn=tmp/(1+eor-tmp)

  # Ghost Probability
  PTitk=fnk*I.Qi.tk
  PTitk=PTitk/Fqn
  PTitk[is.nan(PTitk)]=0
  PTijtk=t(1-apply( cbind(0,PTitk[,-K]) , 1, cumsum))

  # Gradient for a
  wa=Eyn-pn+Emn*(1-pn)
  gradn.a=A*wa

  # Gradient for b

  wb=Eyn*((yn==1)-phaz*Lam0(survn[,2]))+Emn*(1-apply(PTitk*Lnk,1,sum))
  gradn.b=Z*wb

  # Gradient for lambda

  gradn.lam=I.Xieqtk/col.lamk-I.Xi.tk*phaz*Eyn+
    Emn*(PTitk/col.lamk-PTijtk*phaz)


  return(cbind(gradn.a,gradn.b,gradn.lam))
}

# Gradient^2 of log-likelihood
CureCox.Louis.Egradsq=function(survn,yn,A,Z,a,b,Lam0)
{
  tk=knots(Lam0)
  K=length(tk)        # number of unique event times
  nk=table(survn[yn==1,2]) # counting ties, usually all 1
  dn= yn==-1   # censoring indicator
  N=nrow(survn) # Sample size
  cN=sum(dn)   # censor sample size

  # Expectation of gradients
  temp=CureCox.Louis.Egradn(survn,yn,A,Z,a,b,Lam0)
  gradn.a=temp[,1:length(a)]
  gradn.b=temp[,length(a)+1:length(b)]
  gradn.lam=temp[,-(1:(length(a)+length(b)))]

  A=cbind(1,A)

  # Logisitc Probability
  eor=c(exp(A%*%a))
  pn=eor/(eor+1)

  # Hazard Proportion
  phaz=c(exp(Z%*%b))

  # Hazard and Density
  Lnk=phaz%*%t(Lam0(tk))
  Snk=exp(-Lnk)
  lamk=diff(c(0,Lam0(tk)))
  col.lamk=(rep(1,N)%*%t(lamk))

  # Indicators
  I.Qi.tk= (survn[,1]%*%t(rep(1,K))) > (rep(1,N)%*%t(tk))
  I.Xi.tk= (survn[,2]%*%t(rep(1,K))) >= (rep(1,N)%*%t(tk))
  I.Xieqtk= (survn[,2]%*%t(rep(1,K))) == (rep(1,N)%*%t(tk))

  # Expectation on latent variables
  Eyn=yn
  tmp=eor[dn]*exp(-phaz[dn]*Lam0(survn[dn,2]))
  Eyn[dn]=tmp/(1+tmp)

  fnk=outer(c(phaz),lamk)*Snk
  Fqn=apply(fnk*I.Qi.tk,1,sum)
  tmp=eor*Fqn
  Emn=tmp/(1+eor-tmp)

  # Variance on latent variables
  Varyn=rep(0,N)
  Varyn[dn]=Eyn[dn]*(1-Eyn[dn])

  Varmn=Emn*(1+eor)/(1+eor-tmp)
  Emnsq=Varmn+Emn^2

  # Ghost Probability
  PTitk=fnk*I.Qi.tk
  PTitk=PTitk/Fqn
  PTitk[is.nan(PTitk)]=0
  PTijtk=t(1-apply( cbind(0,PTitk[,-K]) , 1, cumsum))

  sumtQPS=-apply(PTitk*Lnk,1,sum)
  sumtkQPS=sumtQPS+t(apply(cbind(0,(PTitk*Lnk)[,-K]),1,cumsum))

  # Grad.aa
  waa=(1-pn)^2*Varmn+Varyn
  grad.aa=t(gradn.a)%*%gradn.a+t(A)%*%diag(waa)%*%A

  # beta terms
  wb.I=(yn==1)-phaz*Lam0(survn[,2])

  # Grad.ab
  wab=Varyn*(wb.I) +
    Varmn*(1-pn)*(1+sumtQPS)
  grad.ab=t(gradn.a)%*%gradn.b + t(A) %*% diag(wab)%*% Z


  # Grad.bb
  wbb=Varyn*(wb.I)^2 +
    Varmn*(1+sumtQPS)^2 +
    Emn*(apply(PTitk*Lnk^2,1,sum)-sumtQPS^2)
  grad.bb=t(gradn.b)%*%gradn.b+t(Z) %*% diag(wbb) %*% Z

  # Lambda terms
  wll.I=I.Xieqtk*(yn==1)/col.lamk- I.Xi.tk*phaz
  wll.II=PTitk/col.lamk-PTijtk*phaz
  wll.III=PTitk/col.lamk^2-2*PTitk/col.lamk*phaz+PTijtk*phaz^2


  # Grad.alam
  walam=wll.I*Varyn+wll.II*Varmn*(1-pn)
  grad.alam=t(gradn.a)%*%gradn.lam+t(A)%*%walam

  # Grad.blam
  wblam= wll.I *  Varyn *(wb.I) +
    wll.II*Varmn*(1+sumtQPS) -
    (PTitk*(Lnk+sumtQPS)/col.lamk-PTijtk*phaz*sumtQPS+sumtkQPS*phaz)*Emn
  grad.blam=t(gradn.b)%*%gradn.lam+t(Z)%*%wblam


  # Grad.lamlam
  wll.offdiag=apply(-(wll.I*Eyn+wll.II*Emn)*phaz,2,sum)
  maxid=pmax(rep(1:K,K),rep(1:K,rep(K,K)))
  wll.block=matrix(wll.offdiag[maxid],K,K)
  diag(wll.block)=apply(wll.I^2*Eyn+wll.III*Emn,2,sum)
  grad.lamlam= t(wll.I)%*%diag(Eyn*Emn)%*%wll.II +
    t(wll.II)%*%diag(Eyn*Emn)%*% wll.I +
    t(wll.II)%*%diag(Emnsq-Emn)%*%wll.II + wll.block

  out=cbind(grad.aa,grad.ab)
  out=rbind(out,cbind(t(grad.ab),grad.bb))
  out=cbind(out,rbind(grad.alam,grad.blam))
  out=rbind(out,cbind(t(grad.alam),t(grad.blam),grad.lamlam))

  return(out)
}

