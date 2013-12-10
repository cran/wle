\name{wle.vonmises}
\alias{wle.vonmises}
\alias{print.wle.vonmises}

\title{von Mises Weighted Likelihood Estimates}

\description{
Computes the weighted likelihood estimates for the parameters of a von Mises distribution:  the mean direction and the concentration parameter.
}

\usage{
wle.vonmises(x, boot = 30, group, num.sol = 1, raf = "HD", smooth, tol =
10^(-6), equal = 10^(-3), max.iter = 500, bias = FALSE, mle.bias =
FALSE, max.kappa = 500, min.kappa = 0.01, use.smooth = TRUE, alpha =
NULL, p = 2, verbose = FALSE, control.circular = list())
\method{print}{wle.vonmises}(x, digits = max(3, getOption("digits") - 3), ...)

}

\arguments{
  \item{x}{a vector. The object is coerced to class \code{\link[circular]{circular}}.}
  \item{boot}{the number of starting points based on boostrap subsamples to use in the search of the roots.}
  \item{group}{the dimension of the bootstap subsamples.}
  \item{num.sol}{maximum number of roots to be searched.}
  \item{raf}{type of Residual adjustment function to be use:

\code{raf="HD"}: Hellinger Distance RAF,

\code{raf="NED"}: Negative Exponential Disparity RAF,

\code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}
  \item{smooth}{the value of the smoothing parameter.}
  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}
  \item{equal}{the absolute value for which two roots are considered the same. (This parameter must be greater than \code{tol}).}
  \item{max.iter}{maximum number of iterations.}
  \item{bias}{logical, if \code{TRUE}, the estimate for kappa is
    computed with a bias corrected method.  Default is \code{FALSE},
    i.e. no bias correction. }
  \item{mle.bias}{logical, if \code{TRUE} a  bias corrected method is
	used to estimate the concentration parameter for the initial values.}
  \item{max.kappa}{maximum value for the concentration parameter.}
  \item{min.kappa}{minimum value for the concentration parameter.}
  \item{use.smooth}{logical, if \code{TRUE} a smoothed model is used,
    default is \code{TRUE}.}
  \item{alpha}{if not \code{NULL} overrides the value of \code{p}. See the next argument \code{p}. This is a different
    parameterization, \code{alpha=-1/2} provides Hellinger Distance RAF,
    \code{alpha=-1} provides Kullback-Leibler RAF and \code{alpha=-2}
    provides Neyman's Chi-Square RAF.}
  \item{p}{this parameter works only when \code{raf="HD"}. \code{p=2}
	provides Hellinger Distance RAF, \code{p=-1}
	provides Kullback-Leibler RAF and \code{p=Inf} provides Neyman's
	Chi-Square RAF.}
  \item{verbose}{logical, if \code{TRUE} warnings are printed.}
  \item{control.circular}{the attribute of the resulting object (\code{mu})}    
  \item{digits}{integer indicating the precision to be used.}
  \item{...}{further parameters in \code{print.wle.vonmises}.}
}

\details{
 Parameters \code{p} and \code{raf} will be change in the future. See
 the reference below for the definition of all the RAF.
}

\value{
  Returns a list with the following components:

  \item{call}{the match.call().}
  \item{mu}{the estimate of the mean direction or the value supplied. If
  \code{num.sol} > 1 then \code{mu} may have length greater than 1, i.e, one value for each root found.}
  \item{kappa}{the estimate of the concentration parameter or the
    value supplied. If \code{num.sol} > 1 then \code{kappa} may have length greater than 1, i.e, one value for each root found.}
  \item{tot.weights}{the sum of the weights divide by the number of
	observations, one value for each root found.}
\item{weights}{the weights associated to each observation, one column vector for each root found.}
\item{f.density}{the non-parametric density estimation.}
\item{m.density}{the smoothed model.}
\item{delta}{the Pearson residuals.}
\item{tot.sol}{the number of solutions found.}
\item{not.conv}{the number of starting points that does not converge after the \code{max.iter} iteration are reached.}
}

\references{
  C. Agostinelli. Robust estimation for circular data.
  Computational Statistics & Data Analysis, 51(12):5867-5875, 2007.
}

\author{Claudio Agostinelli}

\seealso{\code{\link[circular]{circular}}, \code{\link[circular]{mle.vonmises}}.}

\examples{

x <- c(rvonmises(n=50, mu=circular(0), kappa=10), rvonmises(n=5, mu=circular(pi/2), kappa=20))
wle.vonmises(x, smooth=20, group=5)

}

\keyword{robust}

