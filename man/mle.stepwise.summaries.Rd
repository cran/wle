\name{mle.stepwise.summaries}
\title{Accessing summaries for mle.stepwise}
\usage{
summary.mle.stepwise(object, num.max=20, ...)

print.mle.stepwise(x, digits = max(3, getOption("digits") - 3), ...)

print.summary.mle.stepwise(x, digits = max(3, getOption("digits") - 3), ...)
}

\alias{summary.mle.stepwise}
\alias{print.mle.stepwise}
\alias{print.summary.mle.stepwise}

\arguments{
 \item{object}{an object of class \code{mle.stepwise}.}
 \item{x}{an object of class \code{mle.stepwise} or \code{summary.mle.stepwise}.}
 \item{num.max}{the number of the last iterations reported.}}

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





