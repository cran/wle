\name{wle.stepwise.summaries}
\title{Accessing summaries for wle.stepwise}
\usage{
summary.wle.stepwise(object, num.max=20, verbose=FALSE, ...)

print.wle.stepwise(x, digits = max(3, getOption("digits") - 3), ...)

print.summary.wle.stepwise(x, digits = max(3, getOption("digits") - 3), ...)
}

\alias{summary.wle.stepwise}
\alias{print.wle.stepwise}
\alias{print.summary.wle.stepwise}

\arguments{
 \item{object}{an object of class \code{wle.stepwise}.}
 \item{x}{an object of class \code{wle.stepwise} or \code{summary.wle.stepwise}.}
 \item{digits}{number of digits to be used for most numbers.} 
 \item{num.max}{the number of the last iterations reported.}
 \item{verbose}{if \code{TRUE} warnings are printed.}
 \item{...}{additional arguments affecting the summary produced (in
   \code{summary.wle.stepwise}) or further arguments passed to or from
   other methods (in \code{print.wle.stepwise} and
   \code{print.summary.wle.stepwise}).}
}

 
\description{All these functions are \code{\link{methods}} for class \code{wle.stepwise} or \code{summary.wle.stepwise}.}

\details{
The generic accessor functions \code{\link{coefficients}}, \code{\link{fitted.values}}, \code{\link{residuals}} and \code{\link{weights}} can be used to extract various useful features of the value returned by \code{wle.stepwise}. 
}

\value{
  The function \code{summary.wle.stepwise} returns the last \code{num.max} iterations, \code{call} plus:

\item{wstep}{the model for each iteration reported.}
\item{num.max}{the numeber of iterations reported.}
\item{type}{the type of selection procedure used.}
\item{f.in}{the in value.}
\item{f.out}{the out value.}

}

\author{Claudio Agostinelli}

%\examples{}

\keyword{robust}
\keyword{regression}


