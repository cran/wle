\name{wle.onestep}
\title{A One-Step Weighted Likelihood Estimator for Linear model}
\usage{
wle.onestep(ydata, xdata, ini.param, ini.scale, inter=1, raf=1, smooth=0.0320018, num.step=1)
}
%- maybe also `usage' for other functions documented here.
\alias{wle.onestep}
%- Also NEED an `\alias' for EACH other function documented here.
\arguments{
 \item{ydata}{A vector containing the observations from the response variable.}
 \item{xdata}{A matrix containing the observations from the explanatory variables.}
 \item{ini.param}{Starting values for the coefficients.}
 \item{ini.scale}{Starting values for the scale parameters.}
 \item{inter}{If 1, intercept is include in the model.}
 \item{raf}{Type of Residual adjustment function to be use:

\code{raf=1}: Hellinger Distance RAF,

\code{raf=2}: Negative Exponential Disparity RAF,

\code{raf=3}: Chi-Squared Disparity RAF.}
 \item{smooth}{The value of the smoothing parameter.}
 \item{num.step}{Number of the steps.}
}
\description{This function evaluate the One-step weighted likelihood estimator for the regression and scale parameters.}

%\details{}

\value{

\item{coefficients}{The parameters estimator.}
\item{standard.error}{An estimation of the standard error of the parameters estimator.}
\item{scale}{An estimation of the error scale.}
\item{residuals}{The unweighted residuals from the estimated model.}
\item{tot.weights}{The sum of the weights.}
\item{weights}{The weights associated to each observation.}

}

\references{

Agostinelli, C., (1997). A one-step robust estimator based on the weighted likelihood methodology, \emph{ Working Paper n. 1997.16}, Department of Statistics, University of Padova.

Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova.

Agostinelli, C., Markatou, M., (1998). A one-step robust estimator for regression based on the weighted likelihood reweighting scheme, \emph{ Statistics \& Probability Letters}, Vol. 37, n. 4, 341-350.

Agostinelli, C., (1998). Verosimiglianza pesata nel modello di regressione lineare,  \emph{ XXXIX Riunione scientifica della Societ\`a Italiana di Statistica}, Sorrento 1998.

}

\author{Claudio Agostinelli}

%\note{}

\seealso{\link{wle.smooth} an algorithm to choose the smoothing parameter for normal distribution and normal kernel, \link{wle.lm} a function for estimating linear models with normal distribution error and normal kernel.}

\examples{
library(wle)
library(lqs)

data(artificial)

result.lts_lqs(y.artificial~x.artificial, method = "lts")

result.wle_wle.onestep(y.artificial,x.artificial,ini.param=result.lts$coefficients,ini.scale=result.lts$scale[1])
}

\keyword{Weighted Likelihood}