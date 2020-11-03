
Surv.cure=function (time, time2, event, type=c("right","counting"), origin = 0, end = Inf)
{
  if (missing(time))
    stop("Must have a time argument")
  if (!is.numeric(time))
    stop("Time variable is not numeric")
  nn <- length(time)
  ng <- (!missing(time)) + (!missing(time2)) + (!missing(event))

  inputAttributes <- list()
  if (!is.null(attributes(time)))
    inputAttributes$time <- attributes(time)
  if (!missing(time2) && !is.null(attributes(time2)))
    inputAttributes$time2 <- attributes(time2)
  if (!missing(event) && !is.null(attributes(event)))
    inputAttributes$event <- attributes(event)
  if (missing(type) ) {
    if (ng == 1 || ng == 2)
      type <- "right"
    else if (ng == 3)
      type <- "counting"
    else stop("No time variable!")
  }
  else {
    # type <- mtype
    if (ng != 3 && ( type == "counting"))
      stop("Wrong number of args for this type of survival data")
    if (ng != 2 && (type == "right"))
      stop("Wrong number of args for this type of survival data")
  }
  if (ng == 1) {
    if (!is.numeric(time))
      stop("Time variable is not numeric")
    ss <- cbind(time = time - origin, status = sign(end-time))
    type <- "right"
  }
  else if (type == "right" ) {
    if (!is.numeric(time))
      stop("Time variable is not numeric")
    if (missing(event))
      event <- time2
    if (length(event) != nn)
      stop("Time and status are different lengths")

    if (is.logical(event))
      status <- as.numeric(event)
    else if (is.numeric(event)) {
      who2 <- !is.na(event)
      if (max(event[who2]) == 2)
        status <- event - 1
      else status <- event
      temp <- (status == 0 | status == 1)
      status <- ifelse(temp, status, NA)
      if (!all(temp[who2], na.rm = TRUE))
        warning("Invalid status value, converted to NA")
    }
    else stop("Invalid status value, must be logical or numeric")

    status[status == 0] = -1
    status[time >= end] = 0

    ss <- cbind(time = time - origin, status = status)
  }
  else if (type == "counting") {
    if (length(time2) != nn)
      stop("Start and stop are different lengths")
    if (length(event) != nn)
      stop("Start and event are different lengths")
    if (!is.numeric(time))
      stop("Start time is not numeric")
    if (!is.numeric(time2))
      stop("Stop time is not numeric")
    temp <- (time >= time2)
    if (any(temp & !is.na(temp))) {
      time[temp] <- NA
      warning("Stop time must be > start time, NA created")
    }
    temp <- (time >= end)
    if (any(temp & !is.na(temp))) {
      time[temp] <- NA
      warning("Start time must be > end time, NA created")
    }

    if (is.logical(event))
      status <- as.numeric(event)
    else if (is.numeric(event)) {
      who2 <- !is.na(event)
      if (max(event[who2]) == 2)
        status <- event - 1
      else status <- event
      temp <- (status == 0 | status == 1)
      status <- ifelse(temp, status, NA)
      if (!all(temp[who2], na.rm = TRUE))
        warning("Invalid status value, converted to NA")
    }
    else stop("Invalid status value")

    status[status == 0] = -1
    status[time2 >= end] = 0
    ss <- cbind(start = time - origin, stop = time2 - origin,
                status = status)
  }

  dimnames(ss) <- list(NULL, dimnames(ss)[[2]])
  attr(ss, "type") <- type
  if (length(inputAttributes) > 0)
    attr(ss, "inputAttributes") <- inputAttributes
  attr(ss, "origin") <- origin
  attr(ss, "end") <- end
  class(ss) <- "Surv.cure"
  ss
}

print.Surv.cure=function(x, digit = getOption("digits"),...)
{
  object = x

  p.status=ncol(object)

  if(digit >=0 )
    object[,-p.status] = round(object[,-p.status],digits=digit)

  status.char = rep(']',nrow(object))
  status.char[object[,p.status] == 0] = '] cured'
  status.char[object[,p.status] == -1] = '+)'

  out = paste(object[,p.status-1],status.char, sep='')

  if(p.status == 3)
    out = paste('[',object[,1],', ',out, sep='')
  else
    out = gsub(')','',gsub(']', '', out))

  print(out,quote=F)
}


"[.Surv.cure" <- function(x, i, j, ...)
{
  x.mat=x
  class(x.mat)="matrix"

  if(!missing(j))
  {
    return(x.mat[i,j])
  }

  ss = matrix(x.mat[i, ],ncol = ncol(x.mat))
  dimnames(ss) <- dimnames(x)
  attr(ss, "type") <- attr(x, "type")

  class(ss) <- "Surv.cure"
  ss
}
