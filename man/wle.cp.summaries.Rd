\name{wle.cp.summaries}
\title{Summaries and methods for wle.cp}
\usage{
summary.wle.cp(object, num.max=20, ...)

print.wle.cp(x, digits = max(3, getOption("digits") - 3), ...)

print.summary.wle.cp(x, digits = max(3, getOption("digits") - 3), ...)
}

\alias{summary.wle.cp}
\alias{print.wle.cp}
\alias{print.summary.wle.cp}

\arguments{
 \item{object}{an object of class \code{wle.cp}.}
 \item{x}{an object of class \code{wle.cp} or \code{summary.wle.cp}.}
 \item{num.max}{the max number of models should be reported.}
}

\description{All these functions are \code{\link{methods}} for class \code{wle.cp} or \code{summary.wle.cp}.}

\value{
\code{summary.wle.cp} returns a list:

 \item{wcp}{the first \code{num.max} best models with their Weighted Mallows Cp.}
 \item{num.max}{the number of models reported.}
 \item{call}{}
}

\author{Claudio Agostinelli}

\seealso{\code{\link{wle.cp}} a function for evaluate the Weighted Mallows Cp in the linear models.}

%\examples{}

\keyword{Weighted Likelihood}
\keyword{Weighted Mallows Cp}
\keyword{WCp}
\keyword{model selection}

