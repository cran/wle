\name{mle.cp.summaries}
\title{Summaries and methods for mle.cp}
\usage{
\method{summary}{mle.cp}(object, num.max=20, verbose=FALSE, ...)

\method{print}{mle.cp}(x, digits = max(3, getOption("digits") - 3),
num.max=max(1, nrow(x$cp)), ...)

\method{print}{summary.mle.cp}(x, digits = max(3, getOption("digits") - 3), ...)
}

\alias{summary.mle.cp}
\alias{print.mle.cp}
\alias{print.summary.mle.cp}

\arguments{
 \item{object}{an object of class \code{mle.cp}.}
 \item{x}{an object of class \code{mle.cp} or \code{summary.mle.cp}.}
 \item{digits}{number of digits to be used for most numbers.}  
 \item{num.max}{the max number of models should be reported.}
 \item{verbose}{if \code{TRUE} warnings are printed.}
 \item{...}{additional arguments affecting the summary produced (in
   \code{summary.mle.cp}) or further arguments passed to or from other methods
   (in \code{print.mle.cp} and \code{print.summary.mle.cp}).}
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

\keyword{regression}



