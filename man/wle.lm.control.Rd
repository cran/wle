\name{wle.lm.control}
\Rdversion{1.1}
\alias{wle.lm.control}
\title{Auxiliary for Controlling LM Robust Fitting}
\description{
  Auxiliary function as user interface for lm robust fitting.
  %%Typically
  %%only used when calling \code{\link{wle.cv.rb}} (for now, in the future
  %%also \code{\link{wle.lm}} and many others).
}
\usage{
wle.lm.control(nstart = 30, group = NULL, num.sol = 1, 
  raf = c("HD", "NED", "SCHI2", "GKL", "PWD"), tau = 0.1, 
  cutpoint = 0, powerdown = 1, smooth = 0.031, tol = 10^(-6),
  equal = 10^(-3), max.iter = 500, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nstart}{interger. Number of starting points based on boostrap subsamples to use in the search of the roots.}
  \item{group}{integer. Dimension of the bootstap subsamples. The default value is \eqn{max(round(size/2),var+1)} where \eqn{size} is the number of observations and \eqn{var} is the number of predictors.}
  \item{num.sol}{interger. Maximum number of roots to be searched.}
  \item{raf}{type of Residual adjustment function to be used:

\code{raf="HD"}: Hellinger Distance RAF,

\code{raf="NED"}: Negative Exponential Disparity RAF,

\code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF,

\code{raf="GKL"}: Generalized Kullback-Leibler family RAF (see details),

\code{raf="PWD"}: Power Divergence family RAF (see details).}

  \item{tau}{positive real. Used in selecting the member of the RAF family in the case of \code{GKL} or \code{PWD}.}
  \item{cutpoint}{a value in the interval [0,1].}
  \item{powerdown}{a non negative number.}  
  \item{smooth}{the value of the smoothing parameter; used in the evaluation of weights in the case of continuous models.}
  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}
  \item{equal}{the absolute value for which two roots are considered the same. Two roots are compared using the corresponding final weights.}
  \item{max.iter}{maximum number of iterations.}
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
Agostinelli, C., (1998) Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{Ph.D Thesis}, Department of Statistics, University of Padova.

Agostinelli, C., Markatou, M., (1998) A one-step robust estimator for regression based on the weighted likelihood reweighting scheme, \emph{ Statistics \& Probability Letters}, Vol. 37, n. 4, 341-350.

Agostinelli, C., (1998) Verosimiglianza pesata nel modello di regressione lineare,  \emph{ XXXIX Riunione scientifica della Societ\`a Italiana di Statistica}, Sorrento 1998.
}

\author{
Claudio Agostinelli 
}

\keyword{models}
