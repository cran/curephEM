cureph=function (formula, formula2 , data,  subset, na.action, init, control,
                 method = c("EM"), singular.ok = TRUE,
                  var = c("Louis"),...)
{
  method <- match.arg(method)
  var <- match.arg(var)

  Call <- match.call()
  indx1 <- match(c("formula", "data", "subset",
                  "na.action"), names(Call), nomatch = 0)

  if (indx1[1] == 0)
    stop("A formula argument is required")

  temp <- Call[c(1, indx1)]
  temp[[1]] <- as.name("model.frame")
  temp$formula <- if (missing(data))
      terms(formula)
  else
       terms(formula, data = data)
 # temp$formula2 <- if (missing(data))
#    terms(formula2)
#  else terms(formula2, data = data)

  mf1 <- eval(temp, parent.frame())
  if (nrow(mf1) == 0)
    stop("No (non-missing) observations")
  Terms1 <- terms(mf1)


  extraArgs <- list(n.data=nrow(mf1), ...)
  if (length(extraArgs)) {
    controlargs <- names(formals(cureph.control))
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L)
    if (any(indx == 0L))
      stop(gettextf("Argument %s not matched", names(extraArgs)[indx ==
                                                                  0L]), domain = NA)
  }
  if (missing(control))
  {
    control <- cureph.control(n.data=nrow(mf1), ...)
  }


  Y <- model.extract(mf1, "response")
  if (!inherits(Y, "Surv.cure"))
    stop("Response must be a survival-cure object")
  type <- attr(Y, "type")

  if (type != "right" && type != "counting" )
    stop(paste("CurePH model doesn't support \"", type, "\" survival data",
               sep = ""))

    data.n <- nrow(Y)
  if (type == "right")
  {
    time2 = Y[,1]
    event = as.numeric(Y[,2] == 1)
    origin = attributes(Y)$origin
    time = rep(origin,data.n)
    end = min(time2[Y[,2] == 0])

    Y=Surv.cure(time,time2,event, origin=origin, end=end)
  }

  #***** Missing terms.inner function
  #  if (length(attr(Terms, "variables")) > 2) {
  #    ytemp <- terms.inner(attr(Terms, "variables")[1:2])
  #    xtemp <- terms.inner(attr(Terms, "variables")[-2])
  #    if (any(!is.na(match(xtemp, ytemp))))
  #      warning("a variable appears on both the left and right sides of the formula")
  #  }


  contrast.arg <- NULL
  attr(Terms1, "intercept") <- 1
  adrop <- 0

  A <- model.matrix(Terms1, mf1, contrasts = contrast.arg)
  Aatt <- attributes(A)
  xdrop <- Aatt$assign %in% adrop
  A <- A[, !xdrop, drop = FALSE]
  attr(A, "assign") <- Aatt$assign[!xdrop]
  attr(A, "contrasts") <- Aatt$contrasts

  offset1 <- model.offset(mf1)
  if (is.null(offset1) | all(offset1 == 0))
    offset1 <- rep(0, nrow(mf1))
  else if (any(!is.finite(offset1)))
    stop("offsets must be finite")

  assign1 <- survival::attrassign(A, Terms1)
  contr.save1 <- attr(A, "contrasts")

  if(missing(formula2))
  {
    Z=A
    offset2=offset1
    assign2=assign1
    contr.save2=contr.save1
    Terms2=Terms1
  }
  else
  {
    indx2 <- match(c("formula2", "data", "subset",
                     "na.action"), names(Call), nomatch = 0)

    temp <- Call[c(1, indx2)]
    temp[[1]] <- as.name("model.frame")
    temp$formula2 <- if (missing(data))
      terms(formula2)
    else
      terms(formula2, data = data)
    # temp$formula2 <- if (missing(data))
    #    terms(formula2)
    #  else terms(formula2, data = data)
    names(temp)[names(temp)=='formula2']='formula'

    mf2 <- eval(temp, parent.frame())
    Terms2 <- terms(mf2)

    contrast.arg <- NULL
    adrop <- 0

    Z <- model.matrix(Terms2, mf2, contrasts = contrast.arg)
    Zatt <- attributes(Z)
    xdrop <- Zatt$assign %in% adrop
    Z <- Z[, !xdrop, drop = FALSE]
    attr(Z, "assign") <- Zatt$assign[!xdrop]
    attr(Z, "contrasts") <- Zatt$contrasts

    offset2 <- model.offset(mf2)
    if (is.null(offset2) | all(offset2 == 0))
      offset2 <- rep(0, nrow(mf2))
    else if (any(!is.finite(offset2)))
      stop("offsets must be finite")

    assign2 <- survival::attrassign(Z, Terms2)
    contr.save2 <- attr(Z, "contrasts")
  }




  if (missing(init))
    init <- vector('list',3)
  else {
      if(length(init)!= 3)
        stop("wrong length for init argument.
             NPML requires a list: initial vector for alpha, initial vector for beta,
             initial step function for Lambda0. ")
      if(length(init[[1]])!=ncol(A)+1)
        stop("wrong length for init vector for alpha")
    if(length(init[[2]])!=ncol(Z))
      stop("wrong length for init vector for beta")


    temp <- c(cbind(1,A) %*% init[[1]] - sum(c(1,colMeans(A)) * init[[1]]),
              Z %*% init[[2]] - sum(colMeans(Z) * init[[2]])   )
    if (any(temp < .Machine$double.min.exp | temp > .Machine$double.max.exp))
      stop("initial values lead to overflow or underflow of the exp function")
  }

  ###############################################

    #### Add NPML fitter
    if(method == 'EM') {
      fitter <- get("CureCox.EM")
    }
    else if (method == "BFGS-Newton") {
      stop("Not available in this version")
#      fitter <- get("CureCox.NPML")
    }
    else stop(paste("Unknown method", method))

  if(var == 'Louis') {
    var.method <- get("CureCox.Louis.varMatrix")
  }
  else if (var == "Asymptotic") {
    stop("Not available in this version")
    #    var.method <- get("CureCox.asymp.varMatrix")
  }else
      stop(paste("Unknown variance estimation method", method))

    fit <- fitter(Y, A, Z, init, control,var.method)



  if (is.character(fit)) {
    fit <- list(fail = fit)
    class(fit) <- "cureph"
  }
  else {
    if (!is.null(fit$coefficients) && any(is.na(fit$coefficients))) {
      vars <- (1:length(fit$coefficients))[is.na(fit$coefficients)]
      msg <- paste("X matrix deemed to be singular; variable",
                   paste(vars, collapse = " "))
      if (singular.ok)
        warning(msg)
      else stop(msg)
    }
    fit$n <- data.n
    fit$nevent <- sum(Y[, ncol(Y)]==1)
    fit$terms.logistic <- Terms1
    fit$terms.cox <- Terms2
    fit$assign.logistic <- assign1
    fit$assign.cox <- assign2
    class(fit) <- 'cureph'

    if (length(fit$coefficients) && is.null(fit$wald.test)) {
      nabeta <- !is.na(fit$coefficients)
      if (any(sapply(init[1:2],is.null)))
        temp <- fit$coefficients[nabeta]
      else temp <- (fit$coefficients - unlist(init)[1:length(fit$coefficients)])[nabeta]
      fit$wald.test <- survival::coxph.wtest(fit$var[nabeta, nabeta],
                                   temp, control$toler.chol)$test
      fit$df = sum(nabeta)
    }
    na.action <- attr(mf1, "na.action")
    if (length(na.action))
      fit$na.action <- na.action

  }

  if(length(attr(Terms1,"term.labels"))>0)
  {
    fit$var.levels = pmax(sapply(data[,
                                     apply(outer(attr(Terms1,"term.labels"),colnames(data),'==')
                                          ,1,which)]
                                ,nlevels)-1,1)
  }


  fit$formula.logistic <- formula(Terms1)
  fit$formula.cox <- formula(Terms2)
  fit$contrasts.logistic <- contr.save1
  fit$contrasts.cox <- contr.save2
  if (any(offset1 != 0))
    fit$offset.logistic <- offset1
  if (any(offset2 != 0))
    fit$offset.cox <- offset2
  fit$call <- Call
  fit$method <- method

  colnames(fit$var)=rownames(fit$var)=c("logistic: (Intercept)",
                                        rep(NA, ncol(fit$var)-1))
  if(ncol(A)>0)
  {
    binary.A = apply(A, 2,check.binary)
    mean.A = apply(A,2,mean,na.rm=T)
    mean.A[binary.A] = apply(A,2,min,na.rm=T)[binary.A]
  }
  if(ncol(Z)>0)
  {
    binary.Z = apply(Z, 2,check.binary)
    mean.Z = apply(Z,2,mean,na.rm=T)
    mean.Z[binary.Z] = apply(Z,2,min,na.rm=T)[binary.Z]
  }
  if(ncol(A)+ncol(Z) > 0)
  {
    fit$means <- list()
  }
  fitcoef = fit$coefficients
  fit$coefficients=list(logistic = fitcoef[1:(ncol(A)+1)])
  if(ncol(A)>0)
  {
    fit$means$logistic = mean.A
    names(fit$means$logistic) <- colnames(A)
    colnames(fit$var)[1+1:ncol(A)]=rownames(fit$var)[1+1:ncol(A)] =
      paste('logistic:', colnames(A))
  }
  if(ncol(Z)>0)
  {
    fit$means$cox = mean.Z
    fit$coefficients$cox = fitcoef[-(1:(ncol(A)+1))]
    names(fit$coefficients$cox) <- names(fit$means$cox) <- colnames(Z)
    colnames(fit$var)[-1-0:ncol(A)]=rownames(fit$var)[-1-0:ncol(A)] =
      paste('cox:', colnames(Z))
  }
  names(fit$coefficients$logistic) <-  c('(Intercept)', colnames(A))

  fit$data = data
  fit$Y = Y

  fit$origin = attr(Y,'origin')
  fit$end = attr(Y,'end')
  if(is.infinite(fit$end))
  {
    fit$end = max(Y[,2]) + 0.01 * (max(Y[,2])-min(Y[,2]))
  }


  fit
}




