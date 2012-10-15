\name{wle.glm.control}
\Rdversion{1.1}
\alias{wle.glm.control}
\title{Auxiliary for Controlling GLM Robust Fitting}
\description{
Auxiliary function as user interface for glm robust fitting. Typically only used when calling \code{\link{wle.glm}} or \code{\link{wle.glm.fit}}.
}
\usage{
wle.glm.control(boot = 30, group = NULL, num.sol = 1, 
  raf = c("GKL", "PWD", "HD", "NED", "SCHI2"), tau = 0.1, 
  cutpoint = 0, powerdown = 1, delta = NULL, smooth = NULL,
  asy.smooth=0.031, tol = 10^(-6), equal = 10^(-3),
  max.iter = 500, window.size = NULL, use.asymptotic = NULL, 
  use.smooth=TRUE, mle.dispersion = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{boot}{interger. Number of starting points based on boostrap subsamples to use in the search of the roots.}
  \item{group}{integer. Dimension of the bootstap subsamples. The default value is \eqn{max(round(size/2),var+1)} where \eqn{size} is the number of observations and \eqn{var} is the number of predictors.}
  \item{num.sol}{interger. Maximum number of roots to be searched.}
  \item{raf}{type of Residual adjustment function to be used:

\code{raf="GKL"}: Generalized Kullback-Leibler family RAF (see details),

\code{raf="PWD"}: Power Divergence family RAF (see details), 

\code{raf="HD"}: Hellinger Distance RAF,

\code{raf="NED"}: Negative Exponential Disparity RAF,

\code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}

  \item{tau}{positive real. Used in selecting the member of the RAF family in the case of \code{GKL} or \code{PWD}.}
  \item{cutpoint}{a value in the interval [0,1].}
  \item{powerdown}{a non negative number.}  
  \item{delta}{between (0,1). Used in the construction of the weights for the Binomial family.}
  \item{smooth}{the value of the smoothing parameter; used in the evaluation of weights in the case of continuous models.}
  \item{asy.smooth}{the value of the smoothing parameter; used in the evaluation of asymptotic weights. or in the case of continuous models.}  
  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}
  \item{equal}{the absolute value for which two roots are considered the same. Two roots are compared using the corresponding final weights.}
  \item{max.iter}{maximum number of iterations.}
  \item{window.size}{positive real or \code{NULL}. The observations with a distance, in the predictors space, less than this threshold are used to estimate the conditional distribution for a given level of the predictor.}
  \item{use.asymptotic}{interger or \code{NULL}. The minimum number of
  observations for the level of the predictors under which asymptotic
  weights are used.}
  \item{use.smooth}{if \code{TRUE} the smoothed model is used in the
  computation of the Pearson Residuals. For now, the option is used only for
  the Gamma family.}
  \item{mle.dispersion}{if \code{TRUE} the weighted likelihood estimator
  for dispersion is used otherwise the weighted chi-sqaured statistics
  is used.}
  \item{verbose}{if \code{TRUE} warnings are printed.}
}

\details{
The Generalized Kullback-Leibler family RAF is defined as:
\deqn{
\ln(tau*x+1)/tau 
}
for \eqn{tau > 0}.

The Power Divergence family RAF is defined as:
\deqn{
tau*((x + 1)^(1/tau) - 1)
}
for \eqn{0 < tau < Inf} while
\deqn{
\ln(x+1)
}
for \eqn{tau=Inf}.

}
\value{
  A list with the arguments as components.
}

\references{
  Agostinelli, C. and Alqallaf, F. (2009) Robust inference in 
  Generalized Linear Models. Manuscript in preparation.
}

\author{
Claudio Agostinelli and Fatemah Alqallaf
}

\seealso{
\code{\link{wle.glm}}
}

\examples{
### A variation on  example(wle.glm) :

## Annette Dobson's example ...
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
oo <- options(digits = 12) # to see more when tracing :
wle.glm.D93X <- wle.glm(counts ~ outcome + treatment, family=poisson(),
                control=list(glm=glm.control(trace = TRUE), 
                wle=wle.glm.control(raf='GKL', tau=0.15)))
options(oo)
coef(wle.glm.D93X)
}

\keyword{models}
