\name{wle.normal}
\title{Robust Estimation in the Normal Model}
\usage{
wle.normal_function(x, boot=30, group, num.sol=1, raf="HD",
                    smooth=0.003, tol=10^(-6), equal=10^(-3), 
                    max.iter=500)
}
\alias{wle.normal}

\arguments{
  \item{x}{a vector contain the observations.}

  \item{boot}{the number of starting points based on boostrap subsamples to use in the search of the roots.}

  \item{group}{the dimension of the bootstap subsamples. The default value is \eqn{max(round(size/4),2)} where \eqn{size} is the number of observations.}

  \item{num.sol}{maximum number of roots to be searched.}

  \item{raf}{type of Residual adjustment function to be use:

\code{raf="HD"}: Hellinger Distance RAF,

\code{raf="NED"}: Negative Exponential Disparity RAF,

\code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}

  \item{smooth}{the value of the smoothing parameter.}

  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}

  \item{equal}{the absolute value for which two roots are considered the same. (This parameter must be greater than \code{tol}).}

  \item{max.iter}{maximum number of iterations.}
}
\description{
\code{wle.normal} is used to robust estimate the location and the scale parameters via Weighted Likelihood, when the sample is iid from a normal distribution with unknown mean and variance.
}

\value{

\code{wle.normal} returns an object of \code{\link{class}} \code{"wle.normal"}.

  Only print method is implemented for this class.

  The object returned by \code{wle.normal} are:

\item{location}{the estimator of the location parameter, one value for each root found.}
\item{scale}{the estimator of the scale parameter, one value for each root found.}
\item{residuals}{the residuals associated to each observation, one column vector for each root found.}
\item{tot.weights}{the sum of the weights divide by the number of observations, one value for each root found.}
\item{weights}{the weights associated to each observation, one column vector for each root found.}
\item{freq}{the number of starting points converging to the roots.}
\item{call}{the match.call().}
\item{tot.sol}{the number of solutions found.}
\item{not.conv}{the number of starting points that does not converge after the \code{max.iter} iteration are reached.}

}

\references{

Markatou, M., Basu, A. and Lindsay, B.G., (1998). Weighted likelihood estimating equations with a bootstrap root search, \emph{ Journal of the American Statistical Association}, 93, 740-750.

Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova.}

\author{Claudio Agostinelli}

\seealso{\link{wle.smooth} an algorithm to choose the smoothing parameter for normal distribution and normal kernel.}

\examples{
library(wle)

data(cavendish)

result <- wle.normal(cavendish)

result

result <- wle.normal(cavendish,boot=20,num.sol=1)

barplot(result$weights,col=2,xlab="Observations",
        ylab="Weights",ylim=c(0,1),
        names.arg=seq(1:length(result$weights)))

}

\keyword{Weighted Likelihood}
\keyword{robust}
 

