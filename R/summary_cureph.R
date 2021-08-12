summary.cureph <- function(object, combine = T, ...)
{
  aliased <- sapply(coef(object), is.na)
  p <- sum(!unlist(aliased))
  if (p > 0) {

    p.logistic=length(coef(object)[[1]])
    coef.p <- unlist(object$coefficients)
    exp.coef=exp(coef.p)
    exp.negcoef=exp(-coef.p)
    covmat.unscaled <- object$var

    dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))

    var.cf <- diag(covmat.unscaled)
    s.err <- sqrt(var.cf)
    tvalue <- coef.p/s.err

    lower.95 = exp(coef.p-qnorm(.975)*s.err)
    upper.95 = exp(coef.p+qnorm(.975)*s.err)

    dn <- c("coef", 'exp(coef)', "se(coef)")

    pvalue <- 2 * pnorm(-abs(tvalue))
    coef.table <- cbind(coef.p, exp.coef, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), c(dn,
                                                  "z", "Pr(>|z|)"))

    conf.int = cbind(exp.coef, exp.negcoef, lower.95, upper.95)
    dimnames(conf.int) = list(names(coef.p),
                              c('exp(coef)', 'exp(-coef)', 'lower .95', 'upper .95'))

    coef.table= list(logistic=coef.table[1:p.logistic, ,drop=F],
                     cox=coef.table[-(1:p.logistic), ,drop=F ])
    rownames(coef.table$logistic)=names(coef(object)$logistic)
    rownames(coef.table$cox)=names(coef(object)$cox)
    conf.int =list(logistic=conf.int[1:p.logistic,,drop=F  ],
                   cox=conf.int[-(1:p.logistic),,drop=F  ])
    rownames(conf.int$logistic)=names(coef(object)$logistic)
    rownames(conf.int$cox)=names(coef(object)$cox)

    if(length(object$var.levels)==0)
    {
      combine = F
    }

    if(combine & length(object$var.levels)>0)
      if((length(coef(object)$logistic) - 1 != length(coef(object)$cox)) ||
         any(names(coef(object)$logistic)[-1] != names(coef(object)$cox)))
      {
        warning("Combined Wald-test not applicable, as the sets of covariates are different. ")
        combine = F

      }else{
        df.wald=object$var.levels*2
        var.pos=cumsum(object$var.levels)
        var.range.logistic=mapply(':',2+c(0,var.pos[-length(var.pos)]),1+var.pos)
        var.range.cox=sapply(var.range.logistic,'+',var.pos[length(var.pos)])
        var.range = mapply('c',var.range.logistic,var.range.cox)

        chi.wald = unlist(sapply(var.range,function(var.range,var,b)
        {if(any(is.na(b[var.range])))
          NA
          else survival::coxph.wtest(var[var.range,var.range],b[var.range])$test}, covmat.unscaled, coef.p))

        p.wald = 1-pchisq(chi.wald, df.wald)
        comb.table = cbind(chi.wald, df.wald ,p.wald)

        dimnames(comb.table)=list(names(object$var.levels),
                                  c('Wald-chi.square', 'df', 'p-value'))
      }
  }
  else {
    combine = F
    coef.table <- matrix( 0 , 0L, 5L)
    dimnames(coef.table) <- list(NULL, c("coef", 'exp(coef)', "se(coef)",
                                         "z", "Pr(>|z|)"))
    covmat.unscaled <- covmat <- matrix( 0 , 0L, 0L)
  }
  keep <- match(c("call", "terms.logistic","terms.cox", 'n','nevent','df',
                  "contrasts.logistic","contrasts.cox",  "wald.test",
                  "iter", "na.action"), names(object), 0L)
  ans <- c(object[keep], list(coefficients = coef.table,combine=combine,
                              conf.int=conf.int, aliased = aliased,
                              cov = covmat.unscaled))
  if(combine)
    ans$comb.wald=comb.table

  class(ans) <- "summary.cureph"

  ans

}

print.summary.cureph=function(x, digits = max(3, getOption("digits") - 3),
                              signif.stars = getOption("show.signif.stars"),
                              ...)
{

  coef.logistic=x$coefficients$logistic
  coef.cox=x$coefficients$cox



  sig.logistic=rep(' ',nrow(coef.logistic))
  sig.logistic[coef.logistic[,5]<=0.1] = '.'
  sig.logistic[coef.logistic[,5]<=0.05] = '*'
  sig.logistic[coef.logistic[,5]<=0.01] = '**'
  sig.logistic[coef.logistic[,5]<=0.001] = '***'

  coef.logistic[,1:3]=round(coef.logistic[,1:3],digits=digits)
  coef.logistic[,4]=round(coef.logistic[,4],digits=digits-1)
  coef.logistic[sig.logistic!='***',5]=round(coef.logistic[sig.logistic!='***',5],digits=digits-1)
  coef.logistic[sig.logistic=='***',5]=signif(coef.logistic[sig.logistic=='***',5],digits=digits-3)

  if(signif.stars)
  {
    coef.logistic=cbind(coef.logistic,sig.logistic)
    colnames(coef.logistic)[6]=''
  }


  sig.cox=rep(' ',nrow(coef.cox))
  sig.cox[coef.cox[,5]<=0.1] = '.'
  sig.cox[coef.cox[,5]<=0.05] = '*'
  sig.cox[coef.cox[,5]<=0.01] = '**'
  sig.cox[coef.cox[,5]<=0.001] = '***'

  coef.cox[,1:3]=round(coef.cox[,1:3],digits=digits)
  coef.cox[,4]=round(coef.cox[,4],digits=digits-1)
  coef.cox[sig.cox!='***',5]=round(coef.cox[sig.cox!='***',5],digits=digits-1)
  coef.cox[sig.cox=='***',5]=signif(coef.cox[sig.cox=='***',5],digits=digits-3)

  if(signif.stars)
  {
    coef.cox=cbind(coef.cox,sig.cox)
    colnames(coef.cox)[6]=''
  }


  cat('Call:\n')
  print(x$call)
  cat('\n')
  cat('Logistic Model: \n')
  print(formula(x$terms.logistic))
  cat('Cox Model: \n')
  print(formula(x$terms.cox))
  cat('\n')
  cat(paste('  n= ',x$n, ', number of events= ',x$nevent, sep='' ))
  cat('\n \n')
  cat('Logistic:\n')
  print(coef.logistic,quote = FALSE)
  cat('\n')
  if(nrow(coef.cox)>0)
  {
    cat('Cox:\n')
    print(coef.cox,quote = FALSE)
  }
  if(x$combine)
  {
    comb.wald = x$comb.wald

    sig.comb=rep(' ',nrow(comb.wald))
    sig.comb[comb.wald[,3]<=0.1] = '.'
    sig.comb[comb.wald[,3]<=0.05] = '*'
    sig.comb[comb.wald[,3]<=0.01] = '**'
    sig.comb[comb.wald[,3]<=0.001] = '***'

    comb.wald[,1]=round(comb.wald[,1],digits=digits-1)
    comb.wald[sig.comb!='***',3]=round(comb.wald[sig.comb!='***',3],digits=digits-1)
    comb.wald[sig.comb=='***',3]=signif(comb.wald[sig.comb=='***',3],digits=digits-3)

    if(signif.stars)
    {
      comb.wald=cbind(comb.wald,sig.comb)
      colnames(comb.wald)[4]=''
    }

    cat('\n')
    cat('Combined Wald tests: \n')

    print(comb.wald,quote=F)
  }

  cat('--- \n')
  cat(paste('Signif. codes:  0', sQuote('***'), "0.001", sQuote('**'), '0.01', sQuote('*'), '0.05', sQuote('.'), '0.1', sQuote(' '), '1\n'))
  cat('\n')
  cat('Logistic:\n')
  print(x$conf.int$logistic,digits=digits)
  cat('\n')
  if(nrow(coef.cox)>0)
  {
    cat('Cox:\n')
    print(x$conf.int$cox,digits=digits)
    cat('\n')
  }
  cat(paste('Wald test =', signif(x$wald.test,digits=digits), 'on', x$df, 'df, p =',
            signif(1-pchisq(x$wald.test,x$df),digits=digits)))

}
