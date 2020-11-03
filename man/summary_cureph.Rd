\name{summary.cureph}
\alias{summary.cureph}
\alias{print.summary.cureph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing Logistic-CoxPH Cure-Rate Model Fit
}
\description{
These functions are all methods for class \code{cureph} or \code{summary.cureph} objects.
}
\usage{
\method{summary}{cureph}(object,combine = T,...)

\method{print}{summary.cureph}(x, digits = max(3, getOption("digits") - 3),
      signif.stars = getOption("show.signif.stars"),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class "\code{cureph}", usually, a result of a call to \code{\link{cureph}}.
}
  \item{x}{an object of class "\code{summary.glm}", usually, a result of a call to \code{summary.glm}.}
  \item{combine}{logical; if \code{TRUE}, a Wald test table of each variable appearing in both models is returned. }
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If TRUE, "significance stars" are printed for each coefficient.}
  \item{...}{Other parameters to match the generic S3 method.}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cureph}}
}


