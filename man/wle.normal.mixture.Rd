\name{wle.normal.mixture}
\alias{wle.normal.mixture}
\alias{print.wle.normal.mixture}
\title{Robust Estimation in the Normal Mixture Model}

\usage{
wle.normal.mixture(x, m, boot=5, group, num.sol=1, raf="HD", 
                   smooth=0.003, tol=10^(-15), equal=10^(-2), 
                   max.iter=1000, all.comp=TRUE, min.size=0.02,
                   min.weights=0.3, boot.start=10, group.start=3, 
                   tol.start=10^(-6), equal.start=10^(-3),
                   smooth.start=0.003, max.iter.start=100, 
                   max.iter.boot=25, verbose=FALSE)
}

\arguments{
  \item{x}{a vector contain the observations.}
  \item{m}{numbers of components.}
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
  \item{all.comp}{try to find all the components.}
  \item{min.size}{see details}
  \item{min.weights}{see details}
  \item{boot.start}{the number of starting points for the starting process.}
  \item{group.start}{the dimension of the bootstap subsamples in the starting process. The default value is \eqn{max(round(group/4),2)}.}
  \item{tol.start}{the absolute accuracy to be used to achieve convergence of the algorithm in the starting process.}
  \item{equal.start}{the absolute value for which two roots are considered the same in the starting process. (This parameter must be greater than \code{tol.start}).}
  \item{smooth.start}{the value of the smoothing parameter in the starting process.}
  \item{max.iter.start}{maximum number of iterations in the starting process.}
  \item{max.iter.boot}{maximum number of iterations of the starting process.}
  \item{verbose}{if \code{TRUE} warnings are printed.}

}
\description{
\code{wle.normal.mixture} is a preliminary version; it is used to robust estimate the location, scale and proportion parameters via Weighted Likelihood, when the sample is iid from a normal mixture univariate  distribution with known \code{m} number of components.
}

\details{
this function use an iterative procedure to evaluate starting
points. First, using \code{\link{wle.normal}} we try to find the biggest
components, then we discard each observation with weight greater than
\code{min.weights}. The \code{\link{wle.normal}} is run on the remain
observations if the ratio between the number of observations and the
original sample size is greater than \code{min.size}. The convergence of
the algorithm is determined by the difference between two
iterations. This stopping rule could have some problems, as soon as
possible it will replace with the one proposed in Markatou (2000)
pag. 485 (5).
}

\value{

\code{wle.normal.mixture} returns an object of \code{\link{class}} \code{"wle.normal.mixture"}.

  Only print method is implemented for this class.

  The objects returned by \code{wle.normal.mixture} are:

    \item{location}{the estimator of the location parameters, one vector for each root found.}
    \item{scale}{the estimator of the scale parameters, one vector for each root found.}
    \item{pi}{the estimator of the proportion parameters, one vector for each root found.}
    \item{tot.weights}{the sum of the weights, divide by the number of observations, one value for each root found.}
    \item{weights}{the weights associated to each observation, one column vector for each root found.}
    \item{f.density}{the non-parametric density estimation.}
    \item{m.density}{the smoothed model.}
    \item{delta}{the Pearson residuals.}
    \item{freq}{the number of starting points converging to the roots.}
    \item{tot.sol}{the number of solutions found.}
    \item{not.conv}{the number of starting points that does not converge after the \code{max.iter} iteration are reached.}
    \item{call}{the match.call().}

}

\references{ 

Markatou, M., (2000) Mixture models, robustness and the weighted likelihood methodology, \emph{Biometrics}, 56, 483-486.

Markatou, M., (2001) A closer look at the weighted likelihood in the context of mixtures, \emph{Probability and Statistical Models with Applications}, Charalambides, C.A., Koutras, M.V. and Balakrishnan, N. (eds.), Chapman and Hall/CRC, 447-467.
}

\author{Claudio Agostinelli}
%\note{}
%\seealso{}

\examples{
library(wle)
set.seed(1234)
x <- c(rnorm(150,0,1),rnorm(50,15,2))
wle.normal.mixture(x,m=2,group=50,group.start=2,boot=5,num.sol=3)
wle.normal(x,group=2,boot=10,num.sol=3)

}

\keyword{robust}
\keyword{models}
