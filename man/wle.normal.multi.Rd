\name{wle.normal.multi}
\title{Robust Estimation in the Normal Multivariate Model}
\usage{
wle.normal.multi_function(x, boot=30, group, num.sol=1,
                          raf="HD", smooth, tol=10^(-6), 
                          equal=10^(-3), max.iter=500)
}
\alias{wle.normal.multi}
\arguments{
  \item{x}{a matrix contain the observations.}

  \item{boot}{the number of starting points based on boostrap subsamples to use in the search of the roots.}

  \item{group}{the dimension of the bootstap subsamples. The default value is \eqn{max(round(size/4),(var*(var+1)/2+var))} where \eqn{size} is the number of observations and \eqn{var} is the number of variables.}

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
\code{wle.normal.multi} is used to robust estimate the location and the covariance matrix via Weighted Likelihood, when the sample is iid from a normal multivariate distribution with unknown means and variance matrix.
}

\value{

\code{wle.normal.multi} returns an object of \code{\link{class}} \code{"wle.normal.multi"}.

  Only print method is implemented for this class.

  The object returned by \code{wle.normal.multi} are:

\item{location}{the estimator of the location parameters, one vector for each root found.}
\item{variance}{the estimator of the covariance matrix, one matrix for each root found.}
\item{tot.weights}{the sum of the weights divide by the number of observations, one value for each root found.}
\item{weights}{the weights associated to each observation, one column vector for each root found.}
\item{freq}{the number of starting points converging to the roots.}
\item{tot.sol}{the number of solutions found.}
\item{call}{the match.call().}
\item{not.conv}{the number of starting points that does not converge after the \code{max.iter} iteration are reached.}

}

\references{

Markatou, M., Basu, A. and Lindsay, B.G., (1998). Weighted likelihood estimating equations with a bootstrap root search, \emph{ Journal of the American Statistical Association}, 93, 740-750.

Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova.}

\author{Claudio Agostinelli}

\seealso{\link{wle.smooth} an algorithm to choose the smoothing parameter for normal distribution and normal kernel.}

\examples{
library(wle)

data(iris)

smooth <- wle.smooth(dimension=4,costant=4,
                    weight=0.5,interval=c(0.3,0.7))

x.data <- as.matrix(iris[iris[,5]=="virginica",1:4])

result <- wle.normal.multi(x.data,boot=20,group=21,
                           num.sol=3,smooth=smooth$root)

result

result <- wle.normal.multi(x.data,boot=20,group=21,
                           num.sol=1,smooth=smooth$root)

barplot(result$weights,col=2,xlab="Observations",
       ylab="Weights",ylim=c(0,1),
       names.arg=seq(1:length(result$weights)))
}

\keyword{models}
\keyword{multivariate}
\keyword{robust}
 
 