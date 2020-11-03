\name{survpred.cureph}
\alias{survpred.cureph}
\alias{survpred}
\alias{plot.survpred.cureph}
\alias{basehaz.cureph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction for a Logistic-CoxPH Cure-Rate Model
}
\description{
Compute Predictions and Survival Curves from a Logistic-CoxPH Cure-Rate Model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{survpred}{cureph}(object, newdata , time, center = F)
basehaz.cureph(object)

\method{plot}{survpred.cureph}(x, pooled = T, censor = x$censored,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a \code{cureph} object.}
  \item{newdata}{a data frame with the same variable names as those that appear in the \code{cureph} formula. It is also valid to use a vector, if the data frame would consist of a single row.

The curve(s) produced will be representative of a cohort whose covariates correspond to the values in \code{newdata}. Default is the data frame used in the original fit for prediction and  the mean of the covariates for survival curve. See \code{\link{cureph.object}}.}

  \item{time}{a numeric vector; all the time points to be evaluated in the survival table in addition to the observed event times. }

  \item{center}{logical; if \code{TRUE}, the baseline survival is evaluated at the mean of the covariates. Otherwise, it is evaluated at zero. }

  \item{x}{a \code{survpred.cureph} object. }
  \item{pooled}{logical; if \code{TRUE}, the marginal survival curve(s) for the cure-rate model will be plotted. Otherwise, the conditional survival curve(s) in the Cox part will be plotted. }
  \item{censor}{a numeric vector of the observed censoring times. If not \code{NULL}, a cross will mark each censored point on the survival curve. }
  \item{...}{Other parameters to match the generic S3 method.}
}
\details{
The means of factor variables are not computed. Instead, the baseline level is used.

The NPMLE baseline survival is estimated in a Breslow-Nelson-Aalen fashion.
}
\value{
 \item{logistic.linear.predict}{the linear predictions for cure indicator logistic model. }
 \item{logistic.linear.prob}{the estimated probability for NOT being cured. }
 \item{cox.linear.predict}{the linear predictions for event time CoxPH model. }
 \item{cox.cumhaz}{the baseline cummulative hazard for event time CoxPH model. }
 \item{surv.cox}{the marginal survival table for the cure-rate model. }
 \item{surv.cureph}{the conditional survival table for event time CoxPH model. }

 \item{}{The object will also contain the following: \code{censored}, \code{origin} and \code{end}.}
}

\author{
Jue (Marquis) Hou
}

\seealso{
\code{\link{cureph}}
}
\examples{
library(curephEM)

# A simulated data set
data("sim.cureph.data")

# Fit logistic-CoxPH cure-rate model
fit=cureph(Surv.cure(time,time2,event,origin=0,end=20)~Z1+Z2+Z3+Z4,data=sim.cureph.data)

# The survival-prediction object
mysurv = survpred(fit,center=TRUE)

oldpar = par()
par(mfrow=c(1,2))
plot(mysurv)
plot(mysurv,FALSE)
par(oldpar)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
