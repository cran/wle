\name{wle.wrappednormal}
\alias{wle.wrappednormal}
\alias{print.wle.wrappednormal}

\title{Wrapped Normal Weighted Likelihood Estimates}

\description{
Computes the weighted likelihood estimates for the parameters of a
Wrapped Normal distribution:  the mean direction and the concentration
parameter (and the scale parameter).
}

\usage{
wle.wrappednormal(x, mu, rho, sd, K, boot = 30, group, num.sol = 1, raf = "HD",
    smooth = 0.0031, tol = 10^(-6), equal = 10^(-3), min.sd = 0.001,
    min.k = 10, max.iter = 100, use.smooth = TRUE, alpha=NULL, p = 2,
verbose = FALSE, control.circular=list())

\method{print}{wle.wrappednormal}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
	\code{\link[circular]{circular}}. }
  \item{mu}{if a values if provided the parameter is considered known.}
  \item{rho}{if a values if provided the parameter (and \code{sd}) is
	considered known.}
  \item{sd}{if a values if provided the parameter (and \code{rho}) is
	considered known.}
  \item{K}{number of elements used to approximate the density of the
	wrapped normal.}
  \item{boot}{the number of starting points based on boostrap subsamples to use in the search of the roots.}
  \item{group}{the dimension of the bootstap subsamples.}
  \item{num.sol}{maximum number of roots to be searched.}
  \item{raf}{type of Residual adjustment function to be use:

\code{raf="HD"}: Hellinger Distance RAF,

\code{raf="NED"}: Negative Exponential Disparity RAF,

\code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}
  \item{smooth}{the value of the smoothing parameter.}
  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}
  \item{equal}{the absolute value for which two roots are considered the
	same. (This parameter must be greater than \code{tol}).}
  \item{min.sd}{minimum value for the \code{sd} parameter.}
  \item{min.k}{minimum number of elements used to approximate the density of the
	wrapped normal.}
  \item{max.iter}{maximum number of iterations.}
  \item{use.smooth}{logical, if \code{TRUE} a smoothed model is used,
    default is \code{TRUE}.}
  \item{alpha}{see the next argument \code{p}. This is a different
    parameterization, \code{alpha=-1/2} provides Hellinger Distance RAF,
    \code{alpha=-1} provides Kullback-Leibler RAF and \code{alpha=-2}
    provides Neyman's Chi-Square RAF.}
  \item{p}{this parameter works only when \code{raf="HD"}. \code{p=2}
	provide Hellinger Distance RAF, \code{p=-1}
	provide Kullback-Leibler RAF and \code{p=Inf} provide Neyman's
	Chi-Square RAF.}
  \item{verbose}{logical, if \code{TRUE} warnings are printed.}
  \item{control.circular}{the attribute of the resulting objects (\code{mu})}    
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
  \item{rho}{the estimate of the concentration parameter or the
    value supplied. If \code{num.sol} > 1 then \code{rho} may have length
	greater than 1, i.e, one value for each root found.}
  \item{sd}{the estimate of the standard deviation parameter or the
    value supplied. If \code{num.sol} > 1 then \code{sd} may have length
	greater than 1, i.e, one value for each root found.}
    
  \item{tot.weights}{the sum of the weights divide by the number of
	observations, one value for each root found.}
\item{weights}{the weights associated to each observation, one column vector for each root found.}
\item{f.density}{the non-parametric density estimation.}
\item{m.density}{the smoothed model.}
\item{delta}{the Pearson residuals.}
\item{tot.sol}{the number of solutions found.}
\item{not.conv}{the number of starting points that does not converge after the \code{max.iter} iteration are reached.}
}

\references{C. Agostinelli (2003) Robust Estimation for Circular Data, manuscript.}

\author{Claudio Agostinelli}

\seealso{\code{\link[circular]{circular}}, \code{\link[circular]{mle.wrappednormal}}.}

\examples{

if (require(circular)) {
    x <- c(rwrappednormal(n=50, mu=circular(0), sd=1), rwrappednormal(n=5, mu=circular(pi/2), sd=0.5))
    wle.wrappednormal(x, smooth=1/20, group=5)
} else {
   cat("Please, install the package 'circular' in order to use this function.\n")
}

}

\keyword{robust}

