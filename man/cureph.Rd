\name{cureph}
\alias{cureph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NPMLE for Logistic-CoxPH Cure-Rate Model
}
\description{
Fit NPMLE for Logistic-CoxPH Cure-Rate Model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cureph (formula, formula2 , data,  subset, na.action, init, control,
                 method = c("EM"), singular.ok = TRUE,
                  var = c("Louis"),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  a \code{formula} object for cure logistic model,with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{\link{Surv.cure}} function.
}

  \item{formula2}{
  a \code{formula} object for event time CoxPH model
  The Response can be suppressed.
  }

  \item{data}{
  a \code{data.frame} in which to interpret the variables named in the \code{formula} and \code{formula2}, or in the \code{subset} and the \code{weights} argument.
  }

  \item{subset}{
  expression indicating which subset of the rows of data should be used in the fit. All observations are included by default.
  }

  \item{na.action}{
  a missing-data filter function. This is applied to the model.frame after any subset argument has been used. Default is \code{options()$na.action}.
  }

  \item{init}{
  a \code{list} of length 3,  initial values of the iteration. Its compoenents are: a vector of \code{formula} coefficients,a vector of\code{formula2} coefficients and a \code{step.fun} step function of baseline cummulative hazard function. Default initial value is the estimate from naive separate \code{glm} and \code{coxph} models.
  }

  \item{control}{
  Object of class \code{\link{cureph.control}} specifying iteration limit and other control options. Default is \code{cureph.control(...)}.
  }

  \item{method}{
  a character string specifying the method for estimation. The default method is the EM algorithm.
  }

  \item{singular.ok}{
  logical value indicating how to handle collinearity in the model matrix. If \code{TRUE}, the program will automatically skip over columns of the X matrix that are linear combinations of earlier columns. In this case the coefficients for such columns will be \code{NA}, and the variance matrix will contain zeros. For ancillary calculations, such as the linear predictor, the missing coefficients are treated as zeros.
  }

  \item{var}{
  a character string specifying the method for variance estimation. The default method is derived from the Louis formula.
  }

  \item{...}{
  other parameters passed to \code{cureph.control}.
  }
}
\details{
The function can handle left truncation.

\code{"BFGS-Newton"} method optimizes the discretized nonparametric likelihood directly. It has 3 steps. First, it roughly search for a local maximal by BFGS method with default setup. Second, it rescales the parameters proportional to the reciprocal of current gradient in another BFGS optimization. Finally, it searches for a solution to the score equations through a Newton's method with Armijo line search.

The \code{"EM"} method uses a ``ghost copy'' algorithm to accomodate left-truncation.

Both variance estimation methods are asymptotically correct, but they have undercoverage issue in case of a small event number. \code{"Louis"} is recommended among the two.

See \cite{Hou et al. (2016)} for detail.
}

\references{
Hou, J., Chambers, C. and Xu, R. (2016). "A nonparametric maximum likelihood approach for partially
observed cured data with left truncation and right-censoring". To be submitted.

}
\author{
Jue (Marquis) Hou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cureph.object}}, \code{\link{summary.cureph}}, \code{\link{survpred.cureph}}
}
\examples{

# A simulated data set
data('sim.cureph.data')
# Or you may generate another one
sim.cureph.data = cureph.simgen()

# Check out the true coefficients and baseline survival
attributes(sim.cureph.data)

# Fit logistic-CoxPH cure-rate model
fit=cureph(Surv.cure(time,time2,event,origin=0,end=20)~Z1+Z2+Z3+Z4,data=sim.cureph.data)

# Check out the summary
summary(fit)

# Baseline cummulative hazard
basehaz.cureph(fit)

}

