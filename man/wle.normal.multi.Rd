\name{wle.normal.multi}
\title{Robust Estimation in the Normal Multivariate Model}
\usage{
wle.normal.multi_function(data, boot, group, num.sol=1, raf=1, smooth=0.00300996, tol=10^(-6), equal=10^(-3), max.iter=500)
}
%- maybe also `usage' for other functions documented here.
\alias{wle.normal.multi}
%- Also NEED an `\alias' for EACH other function documented here.
\arguments{
 \item{data}{A matrix contain the observations.}
 \item{boot}{The number of starting points based on boostrap subsamples to use in the search of the roots.}
 \item{group}{The dimension of the bootstap subsamples.}
 \item{num.sol}{Maximum number of roots to be searched.}
 \item{raf}{Type of Residual adjustment function to be use:

\code{raf=1}: Hellinger Distance RAF,

\code{raf=2}: Negative Exponential Disparity RAF,

\code{raf=3}: Chi-Squared Disparity RAF.
}
 \item{smooth}{The value of the smoothing parameter.}
 \item{tol}{The absolute accuracy to be used to achieve convergence of the algorithm.
}
 \item{equal}{The absolute value for which two roots are considered the same. (This parameter must be greater than \code{tol}).}
 \item{max.iter}{Maximum number of iterations.}
}
\description{
\code{wle.normal.multi} is used to robust estimate the location and the covariance matrix via Weighted Likelihood Methodology, when the sample is iid from a normal multivariate distribution with unknown means and variance matrix.
}

%\details{}

\value{
\item{location}{The estimator of the location parameters, one vector for each root found.}

\item{variance}{The estimator of the covariance matrix, one matrix for each root found.}

\item{tot.weights}{The sum of the weights, one value for each root found.}

\item{weights}{The weights associated to each observation, one column vector for each root found.}

\item{freq}{The number of starting points converging to the roots.}

\item{tot.sol}{The number of solutions found.}

\item{not.conv}{The number of starting points that does not converge after the \code{max.iter} iteration are reached.}

}

\references{

Markatou, M., Basu, A. and Lindsay, B.G., (1998). Weighted likelihood estimating equations with a bootstrap root search, \emph{ Journal of the American Statistical Association}, 93, 740-750.

Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova.}

\author{Claudio Agostinelli}

%\note{}

\seealso{\link{wle.smooth} an algorithm to choose the smoothing parameter for normal distribution and normal kernel.}

\examples{
library(wle)

data(iris)

smooth_wle.smooth(dimension=4,costant=4,weight=0.5,interval=c(0.3,0.7))

x.data_iris[iris[,5]=="virginica",1:4]

result_wle.normal.multi(x.data,boot=20,group=21,num.sol=3,smooth=smooth$root)

print(result)

barplot(result$weights,col=2,xlab="Observations",ylab="Weights",ylim=c(0,1),names.arg=seq(1:length(result$weights)))
}

\keyword{Weighted Likelihood}
 