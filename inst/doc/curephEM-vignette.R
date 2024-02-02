## ----results= 'hide'----------------------------------------------------------
library(curephEM)
devtools::load_all()
set.seed(531)

## -----------------------------------------------------------------------------
sim.cureph.data = cureph.simgen()

## ----results='markup'---------------------------------------------------------
attr(sim.cureph.data, 'true.coef')

## ----results= 'hide',warning=F------------------------------------------------
fit=cureph(Surv.cure(time,time2,event,origin=0,end=20)~Z1+Z2+Z3+Z4,data=sim.cureph.data)

## ----results= 'hide',warning=F------------------------------------------------
fit2=cureph(Surv.cure(time,time2,event,origin=0,end=20)~Z1+Z2+Z3+Z4,
  formula2 = ~ Z1+Z2,data=sim.cureph.data)

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
mysurv = survpred(fit)

## ----fig.show='hold',fig.cap = "Left: Estimated marginal survival at mean-baselevel. Right: Baseline survival in cox part. "----
plot(mysurv)
plot(mysurv, pooled = F)

