\name{wle.glm.weights}
\Rdversion{1.1}
\alias{wle.glm.weights}
\title{
Weights based on Weighted Likelihood for the GLM model
}
\description{
Evaluate the weights for a given GLM model
}
\usage{
wle.glm.weights(y, x, fitted.values, family = gaussian(),
dispersion = 1, raf = "GKL", tau = 0.1, smooth = NULL,
asy.smooth=0.031, window.size = NULL, use.asymptotic = NULL,
use.smooth=TRUE, tol=10^(-6), dist.method = "euclidean",
cutpoint = 0, powerdown = 1)
}
\arguments{
  \item{y}{
  \code{y} is a vector of observations of length \code{n}.
  }
  \item{x}{
  \code{x} is a design matrix of dimension \code{n * p}.
  }
  \item{fitted.values}{
  the fitted mean values, obtained by transforming the linear predictors 
  by the inverse of the link function. Often obtain as a result of 
  \code{\link{wle.glm.fit} call}.
  }
  \item{family}{a description of the error distribution and link
    function to be used in the model. This can be a character string
    naming a family function, a family function or the result of a call
    to a family function.  (See \code{\link{family}} for details of
    family functions.)}
  \item{dispersion}{value of the dispersion parameter. Used only in the
    Gamma family for now.}
  \item{raf}{type of Residual adjustment function to be used:

    \code{raf="GKL"}: Generalized Kullback-Leibler family RAF (see details),

    \code{raf="PWD"}: Power Divergence family RAF (see details), 

    \code{raf="HD"}: Hellinger Distance RAF,

    \code{raf="NED"}: Negative Exponential Disparity RAF,

    \code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}

  \item{tau}{
  positive real. Used in selecting the member of the RAF family in the 
  case of \code{GKL} or \code{PWD}.
  }
  \item{smooth}{
  the value of the smoothing parameter; used in the case of continuous models.
  }
  \item{asy.smooth}{
  the value of the smoothing parameter; used in the evaluation of 
  asymptotic weights.
  }
  \item{window.size}{positive real or \code{NULL}. The observations with a distance, in the predictors space, less than this threshold are used to estimate the conditional distribution for a given level of the predictor.}
  \item{use.asymptotic}{interger or \code{NULL}. The minimum number of
  observations for the level of the predictors under which asymptotic
  weights are used.}
  \item{use.smooth}{if \code{TRUE} the smoothed model is used in the
  computation of the Pearson Residuals. For now, the option is used only for
  the Gamma family.}
  \item{tol}{the tolerance used in the numerical calculations. For now, the option is used only for the Gamma family.}
  \item{dist.method}{
    distance method passed to \code{\link{dist}} to measure the 
    distance between predictor rows.
  }
  \item{cutpoint}{a value in the interval [0,1].}
  \item{powerdown}{a non negative number.}
}

\value{
A list with two components
  \item{weights}{the weights associated to the observations.}
  \item{asy}{logical. If \code{TRUE} the corresponding weights is evaluated using asymptotic considerations based on Anscombe residuals.}
}

\references{
  Agostinelli, C. and Al-quallaf, F. (2009) Robust inference in 
  Generalized Linear Models. Manuscript in preparation.
}

\author{
  Claudio Agostinelli and Fatemah Al-quallaf
}

\seealso{
  \code{\link{wle.glm}}
}

\examples{
\dontshow{require(utils)
example("wle.glm", echo = FALSE)}
# tau=0.1
wgr.D93 <- extractRoot(wle.glm.D93)
# tau=0.2
w1wgr.D93 <- wle.glm.weights(y = wgr.D93$y, x = wgr.D93$x, fitted.values = wgr.D93$fitted.values, family = wgr.D93$family, raf = "GKL", tau = 0.2, smooth = 0.031, window.size = NULL, use.asymptotic = NULL, dist.method = "euclidean")
# tau=0.3
w2wgr.D93 <- wle.glm.weights(y = wgr.D93$y, x = wgr.D93$x, fitted.values = wgr.D93$fitted.values, family = wgr.D93$family, raf = "GKL", tau = 0.3, smooth = 0.031, window.size = NULL, use.asymptotic = NULL, dist.method = "euclidean")

plot(wgr.D93$wle.weights, ylim=c(0,1), ylab='Weights')
points(w1wgr.D93$weights, col=2)
points(w2wgr.D93$weights, col=3)
legend('bottomright', legend=expression(tau==0.1, tau==0.2, tau==0.3), pch=rep(1,3), col=1:3, inset=0.05)
}

\keyword{robust}
