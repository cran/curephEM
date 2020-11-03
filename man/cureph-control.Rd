\name{cureph.control}
\alias{cureph.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ancillary arguments for controling cureph fits
}
\description{
Assign the ancillary arguments for controling cureph fits
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cureph.control(n.data,eps = 1e-09, toler.chol = .Machine$double.eps^0.75,
                 iter.max = 1000, toler.inf = eps^(1/3),
                 line.search=0.5,init.step = 1/n.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.data}{ Sample size of the data}
  \item{eps}{Iteration continues until the relative change in the log partial likelihood is less than eps. Must be positive.}
  \item{toler.chol}{Tolerance for detection of singularity during a Cholesky decomposion of the variance matrix, i.e., for detecting a redundant predictor variable.}
  \item{iter.max}{Maximum number of iterations to attempt for convergence.}
  \item{toler.inf}{Tolerance criteria for the warning message about a possible infinite coefficient value.}
  \item{line.search}{The scale of shrinkage in each Armijo line search step. }
  \item{init.step}{The initial size of step in Newton's method. }
}
\details{
\code{line.search} and \code{init.step} only matter when \code{method = "BFGS-Newton"}.
}
\value{
a list containing the values of each of the above constants
}
\author{
Jue (Marquis) Hou
}

\seealso{
\code{\link{cureph}}
}
