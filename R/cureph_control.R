cureph.control=function(n.data,
                        eps = 1e-09, toler.chol = .Machine$double.eps^0.75,
                        iter.max = 1000, toler.inf = eps^(1/3),
                        line.search=0.5,init.step = 1/n.data)
{
  if (iter.max < 0)
    stop("Invalid value for iterations")
  if (eps <= 0)
    stop("Invalid convergence criteria")
  if (eps <= toler.chol)
    warning("For numerical accuracy, tolerance should be < eps")
  if (toler.inf <= 0)
    stop("The inf.warn setting must be >0")
  if(line.search <= 0 || line.search >= 1)
    stop("line.search must be > 0 and < 1. ")
  if(init.step <= 0)
    stop("The init.step setting must be >0")
  if(init.step < eps)
    warning("For numerical accuracy, initial step size should be >= eps")

  list(eps = eps, toler.chol = toler.chol, iter.max = as.integer(iter.max),
       toler.inf = toler.inf, line.search=line.search, init.step=init.step)
}
