\name{mle.cv.summaries}
\title{Summaries and methods for mle.cv}
\usage{
\method{summary}{mle.cv}(object, num.max=20, verbose=FALSE, ...)

\method{print}{mle.cv}(x, digits = max(3, getOption("digits") - 3), ...)

\method{print}{summary.mle.cv}(x, digits = max(3, getOption("digits") - 3), ...)
}

\alias{summary.mle.cv}
\alias{print.mle.cv}
\alias{print.summary.mle.cv}

\arguments{
 \item{object}{an object of class \code{mle.cv}.}
 \item{x}{an object of class \code{mle.cv} or \code{summary.mle.cv}.}
 \item{digits}{number of digits to be used for most numbers.}
 \item{num.max}{the max number of models should be reported.}
 \item{verbose}{if \code{TRUE} warnings are printed.}
 \item{...}{additional arguments affecting the summary produced (in
   \code{summary.mle.cv}) or further arguments passed to or from other methods
   (in \code{print.mle.cv} and \code{print.summary.mle.cv}).}

}

\description{All these functions are \code{\link{methods}} for class \code{mle.cv} or \code{summary.mle.cv}.}

\value{
\code{summary.mle.cv} returns a list:

 \item{cv}{the first \code{num.max} best models with their estimated prediction error using CV.}
 \item{num.max}{the number of models reported.}
 \item{call}{}
}

\author{Claudio Agostinelli}

\seealso{\code{\link{mle.cv}} a function for evaluate the Cross-Validation selection criterion for linear models.}

%\examples{}

\keyword{regression}



