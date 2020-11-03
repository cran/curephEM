\name{Surv.cure}
\alias{Surv.cure}
\alias{print.Surv.cure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Survival Object with Cure Portion
}
\description{
Create a Survival Object with Cure Portion
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Surv.cure(time, time2, event, type=c("right","counting"), origin = 0, end = Inf)

\method{print}{Surv.cure}(x, digit = getOption("digits"),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{for right censored data, this is the follow up time. For left truncated data, the first argument is the truncation time.
}
  \item{event}{The status indicator, normally \code{0=alive}, \code{1=dead}. Other choices are \code{TRUE/FALSE (TRUE = death)} or \code{1/2 (2=death)}. Although unusual, the event indicator can be omitted, in which case all subjects are assumed to have an event.}
  \item{time2}{ending time of the counting process data only. Intervals are assumed to be open on the left and closed on the right, \code{(start, end]}. Event indicates whether an event occurred at the end of the interval.}
  \item{type}{character string specifying the type of censoring. Possible values are \code{"right"} and \code{"counting"}.}
  \item{origin}{for counting process data, the hazard function origin. This option was intended to be used in conjunction with a model containing time dependent strata in order to align the subjects properly when they cross over from one strata to another, but it has rarely proven useful.}
  \item{end}{the upper bound for event time. Any subject survives beyond this time is considered as an observed cure. }

  \item{x}{a 'Surv.cure' object. }
  \item{digit}{number of decimal digits in the output. }
  \item{...}{Other parameters to match the generic S3 method.}
}

\seealso{
\code{\link{Surv}}
}
\examples{
# A simulated data set
data("sim.cureph.data")
attach(sim.cureph.data)

# Construct and print a Surv.cure object
Surv.cure(time,time2,event,origin=0,end=20)

# Detach the simulated data
detach(sim.cureph.data)
}
