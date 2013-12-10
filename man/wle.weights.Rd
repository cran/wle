\name{wle.weights}
\alias{wle.weights}
\title{Weights based on Weighted Likelihood for the normal model}
\description{This function evaluated the weights for the vector `x` using the vector `y` in the estimation of the density by the kernel density estimator.}
\usage{
wle.weights(x, y=NULL, smooth=0.0031, sigma2, raf=1, location=FALSE,
  max.iter=1000, tol=10^(-6))
}

\arguments{
  \item{x}{the data set for which the weights would be calculate.}
  \item{y}{the data set used to calculate the weights.}
  \item{smooth}{the value of the smoothing parameter.}
  \item{sigma2}{an estimate of the variance.}
  \item{raf}{type of Residual adjustment function to be use:

\code{raf="HD"}: Hellinger Distance RAF,

\code{raf="NED"}: Negative Exponential Disparity RAF,

\code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}
  \item{location}{if \code{TRUE} the location is estimated. Only available when \code{y=NULL}.}
  \item{max.iter}{maximum number of iterations.}
  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}
}

\value{
  \item{weights}{the weights associated to the \code{x} vector.}
  \item{location}{the location.}
  \item{conv}{\code{TRUE} if the convergence is achived.}
}

\author{Claudio Agostinelli}

\keyword{robust}
