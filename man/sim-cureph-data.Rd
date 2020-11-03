\name{sim.cureph.data}
\alias{sim.cureph.data}
\alias{cureph.simgen}
\docType{data}
\title{
Simulated Data for \code{cureph}
}
\description{
A simulated data for \code{cureph} example.
}
\usage{
data("sim.cureph.data")

cureph.simgen() }
\format{
  A data frame with 200 observations on the following 7 variables.
  \describe{
    \item{\code{time}}{entry time.}
    \item{\code{time2}}{exit time.}
    \item{\code{event}}{binary event indicator.}
    \item{\code{Z1}}{numerical vector.}
    \item{\code{Z2}}{factor with levels '0' and '1'. }
    \item{\code{Z3}}{factor with levels 'A', 'B' and 'C'. }
    \item{\code{Z4}}{factor with levels '0', '1' and '2', colinear with \code{Z2} and \code{Z3}. }
    \item{\code{attr(,"true.coef")}}{the true coefficients. }
    \item{\code{attr(,"true.surv0")}}{the true baseline survival function. }
  }
}
