\name{mle.cp.summaries}
\title{Summaries and methods for mle.cp}
\usage{
summary.mle.cp(object, num.max=20, ...)

print.mle.cp(x, digits = max(3, getOption("digits") - 3), ...)

print.summary.mle.cp(x, digits = max(3, getOption("digits") - 3), ...)
}

\alias{summary.mle.cp}
\alias{print.mle.cp}
\alias{print.summary.mle.cp}

\arguments{
 \item{object}{an object of class \code{mle.cp}.}
 \item{x}{an object of class \code{mle.cp} or \code{summary.mle.cp}.}
 \item{num.max}{the max number of models should be reported.}
}

\description{All these functions are \code{\link{methods}} for class \code{mle.cp} or \code{summary.mle.cp}.}

\value{
\code{summary.mle.cp} returns a list:

 \item{cp}{the first \code{num.max} best models with their Mallows Cp.}
 \item{num.max}{the number of models reported.}
 \item{call}{}
}

\author{Claudio Agostinelli}

\seealso{\code{\link{mle.cp}} a function for evaluate the Mallows Cp.}

%\examples{}

\keyword{Mallows Cp}
\keyword{model selection}

