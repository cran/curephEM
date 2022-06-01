\name{cureph.object}
\alias{cureph.object}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cure-Rate Logistic-CoxPH Regression Object
}
\description{
This class of objects is returned by the \code{cureph} class of functions to represent a fitted cure-rate logistic-CoxPH model. Objects of this class have methods for the functions \code{summary} and \code{survpred}.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coefficients}{the vector of coefficients. If the model is over-determined there will be missing values in the vector corresponding to the redundant columns in the model matrix.
}
  \item{var}{the variance matrix of the coefficients. Rows and columns corresponding to any missing coefficients are set to zero.}
  \item{loglik}{ the log-likelihood with the final values of the coefficients.}
  \item{baseline}{the step-function for baseline cummulative hazard. }
  \item{wald.test}{the Wald test of whether the final coefficients differ from the initial values. }
  \item{iter}{number of iterations used.}
  \item{means}{vector of column means of the X matrix. Subsequent survival curves are adjusted to this value. Factors are set to be the baseline level. }
  \item{n}{the number of observations used in the fit.}
  \item{nevent}{the number of events (usually deaths) used in the fit.}
  \item{method}{the computation method used.}
  \item{na.action}{the na.action attribute, if any, that was returned by the \code{na.action} routine.}

\item{...}{The object will also contain the following: \code{terms.logistic}, \code{terms.cox}, \code{assign.logistic}, \code{assign.cox}, \code{df}, \code{var.levels}, \code{formula.logistic}, \code{formula.cox}, \code{contrasts.logistic}, \code{contrasts.cox}, \code{call}, \code{data}, \code{Y}, \code{origin} and \code{end}.}
}

\seealso{
 \code{\link{cureph}}, \code{\link{summary.cureph}}, \code{\link{survpred.cureph}}.
}
