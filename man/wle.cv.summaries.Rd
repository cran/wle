\name{wle.cv.summaries}
\title{Summaries and methods for wle.cv}
\usage{
summary.wle.cv(object, num.max=20, ...)

print.wle.cv(x, digits = max(3, getOption("digits") - 3), ...)

print.summary.wle.cv(x, digits = max(3, getOption("digits") - 3), ...)
}

\alias{summary.wle.cv}
\alias{print.wle.cv}
\alias{print.summary.wle.cv}

\arguments{
 \item{object}{an object of class \code{wle.cv}.}
 \item{x}{an object of class \code{wle.cv} or \code{summary.wle.cv}.}
 \item{num.max}{the max number of models should be reported.}
}

\description{All these functions are \code{\link{methods}} for class \code{wle.cv} or \code{summary.wle.cv}.}

\value{
\code{summary.wle.cv} returns a list:

 \item{wcv}{the first \code{num.max} best models with their estimated prediction error using WCV.}
 \item{num.max}{the number of models reported.}
 \item{call}{}
}

\author{Claudio Agostinelli}

\seealso{\code{\link{wle.cv}} a function for evaluate the Weighted Cross Validation criterion in the linear models.}

%\examples{}

\keyword{regression}
\keyword{robust}

