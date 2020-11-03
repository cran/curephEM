

fit=cureph(Surv.cure(time,time2,event,origin=0,end=20)~Z1+Z2+Z3+Z4,data=sim.cureph.data)

summary(fit)

mysurv = survpred(fit,center=T)
mysurv$surv.cureph
mysurv$surv.cox

par(mfrow=c(1,2))
plot(mysurv)
plot(mysurv,F)

truelog=c(1,-.63,1,0,0,0)
names(truelog) = names(fit$coefficients$logistic)
truecox=c(-.2,.3,0,0,0)
names(truecox) = names(fit$coefficients$cox)
attr(sim.cureph.data,'true.coef') <- list(
  logistic = truelog,
  cox = truecox
)

attr(sim.cureph.data,'true.surv0') <- function(t) pmin(pmax(1-t,0),1)

save(sim.cureph.data, file="data/sim_cureph_data.Rda")
