\name{binary}
\alias{binary}
\title{Convert decimal base number to binary base}
\description{Convert decimal base number to binary base.}
\usage{
binary(x, dim)
}

\arguments{
  \item{x}{a number in decibal base.}
  \item{dim}{the number of digits, if missing the right number of digits is evaluated.}
}

\value{
  \item{binary}{a vector representing the `x` number in binary base.}
  \item{dicotomy}{the same as `binary` but `TRUE` and `FALSE` instead of 1 and 0.}
}

\note{the elements of `binary` and `dicotomy` are in reverse order.}

\author{Claudio Agostinelli}

\examples{
binary(2)
binary(10,dim=5)
}
\keyword{arith}

