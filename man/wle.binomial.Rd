\name{wle.binomial}
\title{Robust Estimation in the Binomial Model}
\usage{
wle.binomial(x, size, boot=30, group, num.sol=1, raf="HD", 
             tol=10^(-6), equal=10^(-3), max.iter=500)
}
\alias{wle.binomial}
\alias{print.wle.binomial}

\arguments{
  \item{x}{a vector contain the number of success in each \code{size} trials.}

  \item{size}{number of trials.}

  \item{boot}{the number of starting points based on boostrap subsamples to use in the search of the roots.}

  \item{group}{the dimension of the bootstap subsamples. The default value is \eqn{max(round(length(x)/4),2)}.}

  \item{num.sol}{maximum number of roots to be searched.}

  \item{raf}{type of Residual adjustment function to be use:

\code{raf="HD"}: Hellinger Distance RAF,

\code{raf="NED"}: Negative Exponential Disparity RAF,

\code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}

  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}

  \item{equal}{the absolute value for which two roots are considered the same. (This parameter must be greater than \code{tol}).}

  \item{max.iter}{maximum number of iterations.}
}
\description{
\code{wle.binomial} is used to robust estimate the proportion parameters via Weighted Likelihood.
}

\value{

\code{wle.binomial} returns an object of \code{\link{class}} \code{"wle.binomial"}.

  Only print method is implemented for this class.

  The object returned by \code{wle.binomial} are:

\item{p}{the estimator of the proportion parameter, one value for each root found.}
\item{tot.weights}{the sum of the weights divide by the number of observations, one value for each root found.}
\item{weights}{the weights associated to each observation, one column vector for each root found.}
\item{call}{the match.call().}
\item{tot.sol}{the number of solutions found.}
\item{not.conv}{the number of starting points that does not converge after the \code{max.iter} iteration are reached.}

}

\author{Claudio Agostinelli}

%\seealso{}

\examples{
library(wle)

set.seed(1234)

x_rbinom(20,p=0.2,size=10)
wle.binomial(x,size=10)

x_c(rbinom(20,p=0.2,size=10),rbinom(10,p=0.9,size=10))
wle.binomial(x,size=10)

}

\keyword{Weighted Likelihood}
\keyword{robust}
 

