\name{wle.riunif}
\alias{wle.riunif}
\title{The Integer Uniform Distribution}
\description{Generates random deviates from an integer uniform distribution, it is based on \code{\link{runif}}. 
}
\usage{
wle.riunif(x, inf, sup)
}

\arguments{
  \item{x}{number of observations.}
  \item{inf, sup}{lower and upper limits of the distribution.}
}

\author{Claudio Agostinelli}

\seealso{\code{\link{runif}}}

\examples{
wle.riunif(10,0,10)
}
\keyword{distribution}

