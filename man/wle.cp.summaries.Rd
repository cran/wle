\name{wle.cp.summaries}
\title{Summaries and methods for wle.cp}
\usage{
\method{summary}{wle.cp}(object, num.max=20, verbose=FALSE, ...)

\method{print}{wle.cp}(x, digits = max(3, getOption("digits") - 3), ...)

\method{print}{summary.wle.cp}(x, digits = max(3, getOption("digits") - 3), ...)
}

\alias{summary.wle.cp}
\alias{print.wle.cp}
\alias{print.summary.wle.cp}

\arguments{
 \item{object}{an object of class \code{wle.cp}.}
 \item{x}{an object of class \code{wle.cp} or \code{summary.wle.cp}.}
 \item{digits}{number of digits to be used for most numbers.}
 \item{num.max}{the max number of models should be reported.}
 \item{verbose}{if \code{TRUE} warnings are printed.}
 \item{...}{additional arguments affecting the summary produced (in
   \code{summary.wle.cp}) or further arguments passed to or from other methods
   (in \code{print.wle.cp} and \code{print.summary.wle.cp}).}
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

\keyword{regression}
\keyword{robust}

