\name{wle.gamma}
\title{Robust Estimation in the Gamma model}
\usage{
wle.gamma(x, boot=30, group, num.sol=1, raf="HD", smooth=0.008, 
          tol=10^(-6), equal=10^(-3), max.iter=500, 
          shape.int=c(0.01,100), shape.tol=10, 
          use.smooth=TRUE, tol.int)
}
\alias{wle.gamma}
\alias{print.wle.gamma}

\arguments{
  \item{x}{a vector contain the observations.}

  \item{boot}{the number of starting points based on boostrap subsamples to use in the search of the roots.}

  \item{group}{the dimension of the bootstap subsamples. The default value is \eqn{max(round(size/4),var)} where \eqn{size} is the number of observations and \eqn{var} is the number of variables.}

  \item{num.sol}{maximum number of roots to be searched.}

  \item{raf}{type of Residual adjustment function to be use:

\code{raf="HD"}: Hellinger Distance RAF,

\code{raf="NED"}: Negative Exponential Disparity RAF,

\code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}

  \item{smooth}{the value of the smoothing parameter.}

  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}

  \item{equal}{the absolute value for which two roots are considered the same. (This parameter must be greater than \code{tol}).}

  \item{max.iter}{maximum number of iterations.}

  \item{shape.int}{a 2 dimension vector for the interval search of the shape parameter at the first iteration.}
 
  \item{shape.tol}{the length of the interval search used in the successive iterations for the shape parameter.}

  \item{use.smooth}{if \code{FALSE} the unsmoothed model is used. This is usefull when the integration routine does not work well.}

  \item{tol.int}{the absolute accuracy to be used in the integration routine. The default value is \eqn{tol*10^{-4}}.}

}
\description{
\code{wle.gamma} is a preliminary version; it is used to robust estimate the shape and the scale parameters via Weighted Likelihood, when the majority of the data are from a gamma distribution.
}

\details{
The gamma is parametrized as follows (\eqn{\alpha = scale}, \eqn{\omega = shape}): 

             \eqn{f(x) = 1/(\alpha^\omega Gamma(\omega)) x^(\omega-1) e^-(x/\alpha)}

for \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\omega > 0}. 

The function use \code{\link{uniroot}} to solve the estimating equation for \eqn{shape}, hence you can expect error from it. You can use \code{shape.int} and \code{shape.tol} to avoid them. It also use a fortran routine (dqagp) to calculate the smoothed model, i.e., evaluate an integral. Sometime the accuracy is not satisfactory, you can use \code{use.smooth=FALSE} to have an approximate estimation using the model instead of the smoothed model.

The Folded Normal distribution is use as kernel. The bandwith is \eqn{smooth*shape/scale^2}.
}

\value{

\code{wle.gamma} returns an object of \code{\link{class}} \code{"wle.gamma"}.

  Only print method is implemented for this class.

  The object returned by \code{wle.gamma} are:

\item{shape}{the estimator of the shape parameter, one value for each root found.}
\item{scale}{the estimator of the scale parameter, one value for each root found.}
\item{tot.weights}{the sum of the weights divide by the number of observations, one value for each root found.}
\item{weights}{the weights associated to each observation, one column vector for each root found.}
\item{call}{the match.call().}
\item{tot.sol}{the number of solutions found.}
\item{not.conv}{the number of starting points that does not converge after the \code{max.iter} iteration are reached.}

}

\references{

Markatou, M., Basu, A. and Lindsay, B.G., (1998). Weighted likelihood estimating equations with a bootstrap root search, \emph{ Journal of the American Statistical Association}, 93, 740-750.

Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova.}

\author{Claudio Agostinelli}

%\seealso{}

\examples{
library(wle)

set.seed(1234)

x_rgamma(100,2,2)

wle.gamma(x)

x_c(rgamma(30,2,2),rgamma(100,20,20))

wle.gamma(x, boot=10, group=10, num.sol=2)

}

\keyword{models}
\keyword{robust}
 

