\name{wle.fracdiff}
\alias{wle.fracdiff}
%\alias{wle.fracdiff.ao}
%\alias{wle.fracdiff.fitted}
%\alias{wle.fracdiff.sigma2}
%\alias{wle.fracdiff.residuals}
%\alias{wle.fracdiff.equation}
%\alias{wle.fracdiff.solve} 
%\alias{wle.fracdiff.pi.coef}
%\alias{wle.fracdiff.xi.coef}
 
\title{Fit Fractional Models to Time Series - Preliminary Version}
\description{This is a preliminary version of functions for the estimation of the fractional parameter via Weighted Likelihood Estimating Equations and a cassification algorithm. The main function is \code{wle.fracdiff}, the remain functions are for internal use and they should not call by the users. They are not documented here.

}

\usage{wle.fracdiff(x, lower, upper, M, group, na.action=na.fail,
tol=10^(-6), equal=10^(-3), raf="HD", smooth=0.0031, smooth.ao=smooth,
boot=10, num.sol=1, x.init=rep(0,M), use.uniroot=FALSE, max.iter.out=20,
max.iter.in=100, max.iter.step=5000, max.iter.start=max.iter.step,
verbose=FALSE, w.level=0.4, min.weights=0.5, init.values=NULL,
num.max=length(x), include.mean=FALSE, ao.list=NULL, elitist=5,
size.generation=5, size.population=10, type.selection="roulette",
prob.crossover=0.8, prob.mutation=0.02, type.scale="none", scale.c=2)
}

\arguments{
  \item{x}{a univariate time series.}
  \item{lower}{the lower end point of the interval to be searched.}
  \item{upper}{the upper end point of the interval to be searched.}
  \item{M}{the order of the finite memory process used to estimate the d parameter.}
  \item{group}{the dimension of the bootstap subsamples.}
  \item{na.action}{function to be applied to remove missing values.}
  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}
  \item{equal}{the absolute value for which two roots are considered the same. (This parameter must be greater than \code{tol}).}
  \item{raf}{type of Residual adjustment function to be use:
      \code{raf="HD"}: Hellinger Distance RAF,
      \code{raf="NED"}: Negative Exponential Disparity RAF,
      \code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}
  \item{smooth}{the value of the smoothing parameter.}
  \item{smooth.ao}{the value of the smoothing parameter used in the outliers classificaton, default equal to \code{smooth}.}
  \item{boot}{the number of starting points based on boostrap subsamples to use in the search of the roots.}
  \item{num.sol}{maximum number of roots to be searched.}
  \item{x.init}{initial values, a vector with the same length of the \code{M} parameter, or a number, default is 0.}
  \item{use.uniroot}{default: \code{FALSE}, if \code{TRUE} in each step the weighted likelihood estimating equations is solved, otherwise, a maximization is performed on a weighted log-likelihood function with fixed weights. The estimators obtain with the two methods is the same.}
  \item{max.iter.out}{maximum number of iterations in the outer loop.}
  \item{max.iter.in}{maximum number of iterations in the inner loop.}
  \item{max.iter.step}{maximum number of iterations in a step.}
  \item{max.iter.start}{maximum number of iterations in the starting process.}
  \item{verbose}{if \code{TRUE} warnings are printed.}
  \item{w.level}{the threshold used to decide if an observation could be an additive outlier.}
  \item{init.values}{a vector with initial values for the d and the
innovations variance.}
  \item{num.max}{maximum number of observations can be considered as possible additive outliers.}
  \item{include.mean}{Should the model include a mean term? The default is \code{TRUE}.}
  \item{ao.list}{possible list of pattern of additive outliers.}
  \item{min.weights}{see details.}
  \item{size.population}{see details.}
  \item{size.generation}{see details.}
  \item{prob.crossover}{see details.}
  \item{prob.mutation}{see details.}
  \item{type.scale}{see details.}
  \item{type.selection}{see details.}
  \item{elitist}{see details.}
  \item{scale.c}{see details.}
}

\details{
\code{min.weight}: the weighted likelihood equation could have more than one solution. These roots appear for particular situation depending on contamination level and type. We introduce the \code{min.weight} parameter in order to choose only between roots that do not down weight everything. This is not still the optimal solution, and perhaps, in the new release, this part will be change. 

The algorithm used to classify the observations as additive outliers is
a simple genetic algorithm as described in Goldberg (1989). The \code{size.population}, \code{size.generation},
\code{type.selection}, \code{prob.crossover}, \code{prob.mutation}, \code{type.scale}, \code{type.selection}, \code{elitist} and
\code{scale.c} are parameters related to this algorithm. 
}

\value{
   \item{d}{the WLE of the fractional parameter.}
    \item{sigma2}{the WLE of the innovations variance.}
    \item{x.mean}{the WLE of the mean.}
    \item{resid}{the residuals.}
    \item{resid.without.ao}{the residuals with the additive outliers effects.}
    \item{resid.with.ao}{the residuals without the additive outliers effects.}
    \item{x.ao}{the time series without the additive outliers effects.}
    \item{call}{the matched call.}
    \item{weights}{the weights.}
    \item{weights.with.ao}{the weights with the additive outliers effects.}
    \item{weights.without.ao}{the weights without the additive outliers effects}
    \item{tot.sol}{the number of solutions found.}
    \item{not.conv}{the number of starting points that does not converge after the \code{max.iter.out} iteration are reached.}
    \item{ao.position}{the position of the additive outliers.}
}

\references{Agostinelli C., Bisaglia L., (2002) Robust estimation of ARFIMA processes, manuscript.

Goldberg, David E., (1989) Genetic Algorithms in Search, Optimization and Machine Learning. Addison-Wesley Pub. Co. ISBN: 0201157675
}

\author{Claudio Agostinelli}

\examples{
  set.seed(1234)
  resw <- wle.fracdiff(Nile, M=100, include.mean=TRUE, lower=0.01, upper=0.96, group=20)
  resw$d
  resw$sigma2
  resw$x.mean

  x <- Nile
  x[50] <- x[50]+4*sd(x)

  set.seed(1234)
  resw <- wle.fracdiff(x, M=100, include.mean=TRUE, lower=0.01, upper=0.96, group=40)
  resw$d
  resw$sigma2
  resw$x.mean
  resw$ao.position
}


\keyword{robust}
\keyword{ts}
