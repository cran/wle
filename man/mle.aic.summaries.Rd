\name{mle.aic.summaries}
\title{Summaries and methods for mle.aic}
\usage{
summary.mle.aic(object, num.max=20, verbose=FALSE, ...)

print.mle.aic(x, digits = max(3, getOption("digits") - 3), ...)

print.summary.mle.aic(x, digits = max(3, getOption("digits") - 3), ...)
}

\alias{summary.mle.aic}
\alias{print.mle.aic}
\alias{print.summary.mle.aic}

\arguments{
 \item{object}{an object of class \code{mle.aic}.}
 \item{x}{an object of class \code{mle.aic} or \code{summary.mle.aic}.}
 \item{num.max}{the max number of models should be reported.}
 \item{digits}{number of digits to be used for most numbers.} 
 \item{verbose}{if \code{TRUE} warnings are printed.}
 \item{...}{additional arguments affecting the summary produced (in
   \code{summary.mle.aic}) or further arguments passed to or from other methods
   (in \code{print.mle.aic} and \code{print.summary.mle.aic}).}
}

\description{All these functions are \code{\link{methods}} for class \code{mle.aic} or \code{summary.mle.aic}.}

\value{
\code{summary.mle.aic} returns a list:

 \item{aic}{the first \code{num.max} best models with their AIC.}
 \item{num.max}{the number of models reported.}
 \item{call}{}
}

\author{Claudio Agostinelli}

\seealso{\code{\link{mle.aic}} a function for evaluate the Akaike Information Criterion.}

%\examples{}

\keyword{regression}




