\name{wle.onestep}
\title{A One-Step Weighted Likelihood Estimator for Linear model}
\usage{
wle.onestep(formula, data=list(), model=TRUE, x=FALSE, 
            y=FALSE, ini.param, ini.scale, raf="HD", 
            smooth=0.0320018, num.step=1, contrasts=NULL)
}
\alias{wle.onestep}

\arguments{

  \item{formula}{a symbolic description of the model to be fit.
    The details of model specification are given below.}

  \item{data}{an optional data frame containing the variables
    in the model.  By default the variables are taken from
    the environment which \code{wle.stepwise} is called from.}

  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding components of the fit (the model frame, the model matrix, the
    response.)}

  \item{ini.param}{starting values for the coefficients.}

  \item{ini.scale}{starting values for the scale parameters.}

  \item{raf}{type of Residual adjustment function to be use:

\code{raf="HD"}: Hellinger Distance RAF,

\code{raf="NED"}: Negative Exponential Disparity RAF,

\code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}

  \item{smooth}{the value of the smoothing parameter.}

  \item{num.step}{number of the steps.}

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}
}

\description{This function evaluate the One-step weighted likelihood estimator for the regression and scale parameters.}


\value{

\code{wle.onestep} returns an object of \code{\link{class}} \code{"wle.onestep"}.

  Only print method is implemented for this class.

  The object returned by \code{wle.onestep} are:

\item{coefficients}{the parameters estimator.}
\item{standard.error}{an estimation of the standard error of the parameters estimator.}
\item{scale}{an estimation of the error scale.}
\item{residuals}{the unweighted residuals from the estimated model.}
\item{fitted.values}{the fitted values from the estimated model.}
\item{tot.weights}{the sum of the weights divide by the number of observations.}
\item{weights}{the weights associated to each observation.}
\item{call}{the match.call().}
\item{contrasts}{}
\item{xlevels}{}
\item{terms}{the model frame.}
\item{model}{if \code{model=TRUE} a matrix with first column the dependent variable and the remain column the explanatory variables for the full model.}
\item{x}{if \code{x=TRUE} a matrix with the explanatory variables for the full model.}
\item{y}{if \code{y=TRUE} a vector with the dependent variable.}

}

\references{

Agostinelli, C., (1997). A one-step robust estimator based on the weighted likelihood methodology, \emph{ Working Paper n. 1997.16}, Department of Statistics, University of Padova.

Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova.

Agostinelli, C., Markatou, M., (1998). A one-step robust estimator for regression based on the weighted likelihood reweighting scheme, \emph{ Statistics \& Probability Letters}, Vol. 37, n. 4, 341-350.

Agostinelli, C., (1998). Verosimiglianza pesata nel modello di regressione lineare,  \emph{ XXXIX Riunione scientifica della Societ\`a Italiana di Statistica}, Sorrento 1998.

}

\author{Claudio Agostinelli}

\seealso{\link{wle.smooth} an algorithm to choose the smoothing parameter for normal distribution and normal kernel, \link{wle.lm} a function for estimating linear models with normal distribution error and normal kernel.}

\examples{
library(wle)
library(lqs)

data(artificial)

result.lts <- lqs(y.artificial~x.artificial, 
                 method = "lts")

result.wle <- wle.onestep(y.artificial~x.artificial,
               ini.param=result.lts$coefficients,
               ini.scale=result.lts$scale[1])

result.wle
}


\keyword{robust}
\keyword{regression}

