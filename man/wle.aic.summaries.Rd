\name{wle.aic.summaries}
\title{Summaries and methods for wle.aic}
\usage{
summary.wle.aic(object, num.max=20, ...)

print.wle.aic(x, digits = max(3, getOption("digits") - 3), ...)

print.summary.wle.aic(x, digits = max(3, getOption("digits") - 3), ...)
}

\alias{summary.wle.aic}
\alias{print.wle.aic}
\alias{print.summary.wle.aic}

\arguments{
 \item{object}{an object of class \code{wle.aic}.}
 \item{x}{an object of class \code{wle.aic} or \code{summary.wle.aic}.}
 \item{num.max}{the max number of models should be reported.}
}

\description{All these functions are \code{\link{methods}} for class \code{wle.aic} or \code{summary.wle.aic}.}

\value{
\code{summary.wle.aic} returns a list:

 \item{waic}{the first \code{num.max} best models with their Weighted Akaike Information Criterion.}
 \item{num.max}{the number of models reported.}
 \item{call}{}
}

\author{Claudio Agostinelli}

\seealso{\code{\link{wle.aic}} a function for evaluate the Weighted Akaike Information Criterion in the linear models.}

%\examples{}

\keyword{Weighted Likelihood}
\keyword{Weighted Akaike Information Criterion}
\keyword{WAIC}
\keyword{model selection}

