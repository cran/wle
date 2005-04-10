\name{wle.aic.ar.summaries}
\title{Summaries and methods for wle.aic.ar}
\usage{
\method{summary}{wle.aic.ar}(object, num.max=20, verbose=FALSE, ...)

\method{print}{wle.aic.ar}(x, digits = max(3, getOption("digits") - 3),
num.max=max(1, nrow(x$waic)), ...)

\method{print}{summary.wle.aic.ar}(x, digits = max(3, getOption("digits") - 3), ...)
}

\alias{summary.wle.aic.ar}
\alias{print.wle.aic.ar}
\alias{print.summary.wle.aic.ar}

\arguments{
 \item{object}{an object of class \code{wle.aic.ar}.}
 \item{x}{an object of class \code{wle.aic} or \code{summary.wle.aic.ar}.}
 \item{digits}{number of digits to be used for most numbers.}
 \item{num.max}{the max number of models should be reported.}
 \item{verbose}{if \code{TRUE} warnings are printed.}
 \item{...}{additional arguments affecting the summary produced (in
   \code{summary.wle.aic.ar}) or further arguments passed to or from other methods
   (in \code{print.wle.aic.ar} and \code{print.summary.wle.aic.ar}).}
}

\description{All these functions are \code{\link{methods}} for class \code{wle.aic.ar} or \code{summary.wle.aic.ar}.}

\value{
\code{summary.wle.aic.ar} returns a list:

 \item{waic}{the first \code{num.max} best models with their Weighted Akaike Information Criterion.}
 \item{num.max}{the number of models reported.}
 \item{call}{}
}

\author{Claudio Agostinelli}

\seealso{\code{\link{wle.aic.ar}} a function for evaluate the Weighted
  Akaike Information Criterion for autoregressive models.}

\keyword{ts}
\keyword{robust}
