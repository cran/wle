\name{wle.poisson}
\title{Robust Estimation in the Poisson Model}
\usage{
wle.poisson(x, boot=30, group, num.sol=1, raf="HD", 
            tol=10^(-6), equal=10^(-3), 
            max.iter=500, verbose=FALSE)
}
\alias{wle.poisson}
\alias{print.wle.poisson}

\arguments{
  \item{x}{a vector contain the number of success.}

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

  \item{verbose}{if \code{TRUE} warnings are printed.}
}
\description{
\code{wle.poisson} is used to robust estimate the lambda parameters in the poisson model via Weighted Likelihood.
}

\value{

\code{wle.poisson} returns an object of \code{\link{class}} \code{"wle.poisson"}.

  Only print method is implemented for this class.

  The object returned by \code{wle.poisson} are:

\item{lambda}{the estimator of the lambda parameter, one value for each root found.}
\item{tot.weights}{the sum of the weights divide by the number of observations, one value for each root found.}
\item{weights}{the weights associated to each observation, one column vector for each root found.}
\item{f.density}{the non-parametric density estimation.}
\item{m.density}{the smoothed model.}
\item{delta}{the Pearson residuals.}
\item{call}{the match.call().}
\item{tot.sol}{the number of solutions found.}
\item{not.conv}{the number of starting points that does not converge after the \code{max.iter} iteration are reached.}

}

\author{Claudio Agostinelli}

\references{

Markatou, M., Basu, A., and Lindsay, B.G., (1997) Weighted likelihood estimating equations: The discrete case with applications to logistic regression, \emph{Journal of Statistical Planning and Inference}, 57, 215-232.

Agostinelli, C., (1998) Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova.}


}

%\seealso{}

\examples{
library(wle)

set.seed(1234)

x <- rpois(40,5)
wle.poisson(x)

x <- c(rpois(40,5),rpois(10,20))
wle.poisson(x)

}

\keyword{models}
\keyword{robust}
 

