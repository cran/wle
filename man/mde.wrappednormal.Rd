\name{mde.wrappednormal}
\alias{mde.wrappednormal}
\alias{print.mde.wrappednormal}

\title{Wrapped Normal Minimum Distance Estimates}
\description{
Computes the minimum distance estimates for the parameters of a Wrapped Normal distribution: the mean direction and the concentration
parameter (and the scale parameter).
}
\usage{
mde.wrappednormal(x, bw, mu = NULL, rho = NULL, sd = NULL, alpha = NULL, p = 2, tol = 1e-05, n = 512, from = circular(0), to = circular(2 * pi), lower = NULL, upper = NULL, method = "L-BFGS-B", lower.rho = 1e-06, upper.rho = 1 - 1e-06, min.sd = 0.001, K = NULL, min.k = 10, control.circular = list(), ...)
\method{print}{mde.wrappednormal}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{a vector. The object is coerced to class \code{\link[circular]{circular}}.}
  \item{bw}{the value of the smoothing parameter.}
  \item{mu}{initial value for the mean direction. Default: maximum likelihood estimate.}
  \item{rho}{initial value for the concentration parameter. Default: maximum likelihood estimate.}
  \item{sd}{initial value for the standard deviation parameter. This value is used only if \code{rho} is \code{NULL}. Default: maximum likelihood estimate.}
   \item{alpha}{see the next argument \code{p}. This is a different
    parameterization, \code{alpha=-1/2} provides Hellinger distance,
    \code{alpha=-1} provides Kullback-Leibler distance and \code{alpha=-2}
    provides Neyman's Chi-Square distance.}
  \item{p}{\code{p=2} provides Hellinger distance, \code{p=-1}
	provides Kullback-Leibler distance and \code{p=Inf} provides Neyman's
	Chi-Square distance. It is ignored if \code{alpha} is not \code{NULL}.}
  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm. This argument is passed to the function which determined the Maximum Likelihood estimates of the parameters. See \code{\link[circular]{mle.wrappednormal}}.}
  \item{n}{number of points used to approximate the density.}
  \item{from}{from which point in the circle the density is approximate.}
  \item{to}{to which point in the circle the density is approximate.}
  \item{lower}{a 2 elements vector passed to \code{optim} used to constrained optimization. First element for the mean direction, second element for the concentration.}
  \item{upper}{a 2 elements vector passed to \code{optim} used to constrained optimization. First element for the mean direction, second element for the concentration.}
  \item{method}{passed to \code{optim}.}
  \item{lower.rho}{if \code{lower} is \code{NULL} this parameter is used to constrained optimization for the concentration parameter.}
  \item{upper.rho}{if \code{upper} is \code{NULL} this parameter is used to constrained optimization for the concentration parameter.}
  \item{min.sd}{minimum value for the \code{sd} parameter. This argument is passed to the function which determined the Maximum Likelihood estimates of the parameters. See \code{\link[circular]{mle.wrappednormal}}.}
  \item{K}{number of elements used to approximate the density of the wrapped normal.}
  \item{min.k}{minimum number of elements used to approximate the density of the wrapped normal.}
  \item{control.circular}{the attribute of the resulting object (\code{mu})}
  \item{digits}{integer indicating the precision to be used.}
  \item{\dots}{further parameters in \code{print.mde.wrappednormal}.}
}
\details{
The distance from an estimated density (by the non parametric kernel density estimator) and the model is evaluated by simple rectangular approximation. \code{optim} is used to performs minimization.

}
\value{
  Returns a list with the following components:

  \item{call}{the match.call().}
  \item{mu}{the estimate of the mean direction.}
  \item{rho}{the estimate of the concentration parameter.}
  \item{sd}{the estimate of the standard deviation parameter.}
  \item{dist}{the distance between the estimated density and the model.}
  \item{data}{the original supplied data converted in radians, clockwise and zero at 0.}
  \item{x}{the 'n' coordinates of the points where the density is estimated.}
  \item{y}{the estimated density values.}
  \item{k}{the density at the model.}
}

\references{C. Agostinelli (2006) Robust Estimation for Circular Data, under revision.}

\author{Claudio Agostinelli}

\seealso{\code{\link[circular]{circular}}, \code{\link[circular]{mle.wrappednormal}} and \code{\link[wle]{wle.wrappednormal}}.}

\examples{

if (require(circular)) {
    set.seed(1234)
    x <- c(rwrappednormal(n=200, mu=circular(0), sd=0.6), rwrappednormal(n=20, mu=circular(pi/2), sd=0.1))
    res <- mde.wrappednormal(x, bw=0.08, mu=circular(0), sd=0.6)
    res
    plot(circular(0), type='n', xlim=c(-1, 1.75), shrink=1.2)
    lines(circular(res$x), res$y)
    lines(circular(res$x), res$k, col=2)
    legend(1,1.5, legend=c('estimated density', 'MDE'), lty=c(1, 1), col=c(1, 2))
} else {
   cat("Please, install the package 'circular' in order to use this function.\n")
}

}

\keyword{robust}

