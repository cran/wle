\name{wle.aic.summaries}
\title{Summaries and methods for wle.aic}
\usage{
summary.wle.aic(object, num.max=20, verbose=FALSE, ...)

print.wle.aic(x, digits = max(3, getOption("digits") - 3), ...)

print.summary.wle.aic(x, digits = max(3, getOption("digits") - 3), ...)
}

\alias{summary.wle.aic}
\alias{print.wle.aic}
\alias{print.summary.wle.aic}

\arguments{
 \item{object}{an object of class \code{wle.aic}.}
 \item{x}{an object of class \code{wle.aic} or \code{summary.wle.aic}.}
 \item{digits}{number of digits to be used for most numbers.}
 \item{num.max}{the max number of models should be reported.}
 \item{verbose}{if \code{TRUE} warnings are printed.}
 \item{...}{additional arguments affecting the summary produced (in
   \code{summary.wle.aic}) or further arguments passed to or from other methods
   (in \code{print.wle.aic} and \code{print.summary.wle.aic}).}
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

\keyword{regression}
\keyword{robust}
