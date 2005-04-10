\name{wle.aic.ar}
\alias{wle.aic.ar}
\alias{wle.ar.wls}
\title{Weighted Akaike Information Criterion for AR models}
\description{
The function evaluate the Weighted Akaike Information Criterion for
AutoRegressive Models. This is a robust model selection method to choose
the order of an AutoRegressive model.
}
\usage{
wle.aic.ar(x, order = c(1, 0), seasonal = list(order = c(0, 0), period = NA), group, group.start, group.step = group.start, xreg = NULL, include.mean = TRUE, na.action = na.fail, tol = 10^(-6), tol.step = tol, equal = 10^(-3), equal.step = equal, raf = "HD", var.full = 0, smooth = 0.0031, smooth.ao = smooth, boot = 10, boot.start = 10, boot.step = boot.start, num.sol = 1, x.init = 0, x.seasonal.init = 0, max.iter.out = 20, max.iter.in = 50, max.iter.start = 200, max.iter.step = 500, verbose = FALSE, w.level = 0.4, min.weights = 0.5, population.size = 10, population.choose = 5, elements.random = 2, wle.start = FALSE, init.values = NULL, num.max = NULL, num.sol.step = 2, min.weights.aic = 0.5, approx.w = TRUE, ask = TRUE, alpha = 2, method = "WLS")
}

\arguments{
  \item{x}{a univariate time series.}
  \item{order}{maximum order to investigate. A specification of the non-seasonal part of the ARI model: the two components \code{(p,d)} are the AR order and the degree of differencing.}
  \item{seasonal}{a specification of the seasonal part of the ARI model, plus the period (which defaults to \code{frequency(x)}).}
  \item{group}{the dimension of the bootstap subsamples.}
  \item{group.start}{the dimension of the bootstap subsamples used in the starting process if \code{wle.init=TRUE}.}
  \item{group.step}{the dimension of the bootstap subsamples used in a step, it must be less than \code{group}.}
  \item{xreg}{optionally, a vector or matrix of external regressors, which must have the same number of rows as \code{x}.}
  \item{include.mean}{Should the ARI model include a mean term? The
default is \code{TRUE} for undifferenced series, \code{FALSE} for differenced ones (where a mean would not affect the fit nor predictions).}
  \item{na.action}{function to be applied to remove missing values.}
  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}
  \item{tol.step}{the absolute accuracy to be used to achieve convergence in a step.}
  \item{equal}{the absolute value for which two roots are considered the same. (This parameter must be greater than \code{tol}).}
  \item{equal.step}{the absolute value for which two roots are considered the same in a step. (This parameter must be greater than \code{tol.step}).}
  \item{raf}{type of Residual adjustment function to be use:
      \code{raf="HD"}: Hellinger Distance RAF,
      \code{raf="NED"}: Negative Exponential Disparity RAF,
      \code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}
  \item{var.full}{An estimate of the residual variance for the full model.}
  \item{smooth}{the value of the smoothing parameter.}
  \item{smooth.ao}{the value of the smoothing parameter used in the outliers classificaton, default equal to \code{smooth}.}
  \item{boot}{the number of starting points based on boostrap subsamples to use in the search of the roots.}
  \item{boot.start}{the number of starting points based on boostrap subsamples to use in the search of the roots in the starting process.}
  \item{boot.step}{the number of starting points based on boostrap subsamples to use in the search of the roots in a step.}
  \item{num.sol}{maximum number of roots to be searched.}
  \item{x.init}{initial values, a vector with the same length of the AR order, or a number, default is 0.}
  \item{x.seasonal.init}{initial values, a vector with the same length of the SAR order, or a number, default is 0.}
  \item{max.iter.out}{maximum number of iterations in the outer loop.}
  \item{max.iter.in}{maximum number of iterations in the inner loop.}
  \item{max.iter.start}{maximum number of iterations in the starting process.}
  \item{max.iter.step}{maximum number of iterations in a step.}
  \item{verbose}{if \code{TRUE} warnings are printed.}
  \item{w.level}{the threshold used to decide if an observation could be an additive outlier.}
  \item{min.weights}{see details.}
  \item{population.size}{see details.}
  \item{population.choose}{see details.}
  \item{elements.random}{see details.}
  \item{wle.start}{if \code{TRUE} a weighted likelihood estimation is used to have a starting value.}
  \item{init.values}{a vector with initial values for the AR and seasonal AR coefficients and the innovations variance.}
  \item{num.max}{maximum number of observations can be considered as possible additive outliers.}
  \item{num.sol.step}{maximum number of roots to be searched in a step.}
  \item{min.weights.aic}{see details.}
  \item{approx.w}{logical: if \code{TRUE} an approximation is used to
    evaluate the weights in the outlier identification procedure.}
  \item{ask}{logical. If \code{TRUE}, in the case of multiple roots in
    the full model, the users is asked for selecting the root.}
  \item{alpha}{penalty value.}
  \item{method}{if "WLE" the parameters are estimated using weighted
    likelihood estimating equations in the reduced models, otherwise if
    "WLS" a weighted least squares approach is used with weights based
    on the full model.}
}
\details{
\code{min.weights}: the weighted likelihood equation could have more than one solution. These roots appear for particular situation depending on contamination level and type. We introduce the \code{min.weight} parameter in order to choose only between roots that do not down weight everything. This is not still the optimal solution, and perhaps, in the new release, this part will be change. 

\code{min.weights.aic} is used as \code{min.weights} but in the full model.
The algorithm used to classify the observations as additive outliers is
made by a genetic algorithm. The \code{population.size},
\code{population.choose} and \code{elements.random} are parameters
related to this algorithm.

The function \code{wle.ar.wls} is used to estimate the parameter of an
autoregressive model by weighted least squares where the weights are
those from the weighted likelihood estimating equation of the full
model (the model with the hightest order). 

}

\value{

  A list of class \code{wle.aic.ar} with the following components:
  \item{full.model}{the results for the full model, that is an object of class \code{wle.arima} see \code{\link{wle.ar}} help for further details.} 
  \item{waic}{Weighted Akaike Information Criterion for each submodels.}
  \item{call}{\code{match.call} result.}
}

\references{

 Agostinelli C, (2004) Robust Akaike Information Criterion
 for ARMA models, Rendiconti per gli Studi Economici Quantitativi,
 1-14, isbn: 88-88037-10-1.

 Agostinelli C., (2003) Robust time series estimation via weighted
likelihood, in: Development in Robust Statistics. International
Conference on Robust Statistics 2001, Eds. Dutter, R. and Filzmoser,
P. and Rousseeuw, P. and Gather, U., Physica Verlag.

}

\author{Claudio Agostinelli}

\seealso{\code{\link{wle.ar}}}
\examples{

data(rocky)

res <- wle.aic.ar(x=rocky, order=c(6,0), group=50, group.start=30, method="WLS")
res
plot(res$full.model$weights)
}

\keyword{robust}
\keyword{ts}
