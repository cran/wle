\name{mle.stepwise.summaries}
\title{Accessing summaries for mle.stepwise}
\usage{
\method{summary}{mle.stepwise}(object, num.max=20, verbose=FALSE, ...)

\method{print}{mle.stepwise}(x, digits = max(3, getOption("digits") - 3), max(1, nrow(x$step))...)

\method{print}{summary.mle.stepwise}(x, digits = max(3, getOption("digits") - 3), ...)
}

\alias{summary.mle.stepwise}
\alias{print.mle.stepwise}
\alias{print.summary.mle.stepwise}

\arguments{
 \item{object}{an object of class \code{mle.stepwise}.}
 \item{x}{an object of class \code{mle.stepwise} or
   \code{summary.mle.stepwise}.}
 \item{digits}{number of digits to be used for most numbers.}
 \item{num.max}{the number of the last iterations reported.}
 \item{verbose}{if \code{TRUE} warnings are printed.}
 \item{...}{additional arguments affecting the summary produced (in
   \code{summary.mle.stepwise}) or further arguments passed to or from other methods
   (in \code{print.mle.stepwise} and \code{print.summary.mle.stepwise}).} 
}

\description{All these functions are \code{\link{methods}} for class \code{mle.stepwise} or \code{summary.mle.stepwise}.}

%\details{}

\value{
  The function \code{summary.mle.stepwise} returns the last \code{num.max} iterations, \code{call} plus:

\item{step}{the model for each iteration reported.}
\item{num.max}{the number of iterations reported.}
\item{type}{the type of selection procedure used.}
\item{f.in}{the in value}
\item{f.out}{the out value}

}

\author{Claudio Agostinelli}

%\examples{}

\keyword{regression}





