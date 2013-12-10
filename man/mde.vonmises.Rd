\name{mde.vonmises}
\alias{mde.vonmises}
\alias{print.mde.vonmises}

\title{von Mises Minimum Distance Estimates}
\description{
Computes the minimum distance estimates for the parameters of a von Mises distribution:  the mean direction and the concentration parameter.
}
\usage{
mde.vonmises(x, bw, mu = NULL, kappa = NULL, n = 512,
  from = circular(0), to = circular(2 * pi), lower = NULL,
  upper = NULL, method = "L-BFGS-B", lower.kappa = .Machine$double.eps,
  upper.kappa = Inf, alpha = NULL, p = 2, control.circular = list(), ...)
\method{print}{mde.vonmises}(x, digits = max(3, getOption("digits") - 3), ...)

}

\arguments{
  \item{x}{a vector. The object is coerced to class \code{\link[circular]{circular}}.}
  \item{bw}{the value of the smoothing parameter.}
  \item{mu}{initial value for the mean direction. Default: maximum likelihood estimate.}
  \item{kappa}{initial value for the concentration parameter. Default: maximum likelihood estimate.}
  \item{n}{number of points used to approximate the density.}
  \item{from}{from which point in the circle the density is approximate.}
  \item{to}{to which point in the circle the density is approximate.}
  \item{lower}{a 2 elements vector passed to \code{optim} used to constrained optimization. First element for the mean direction, second element for the concentration.}
  \item{upper}{a 2 elements vector passed to \code{optim} used to constrained optimization. First element for the mean direction, second element for the concentration.}
  \item{method}{passed to \code{optim}.}
  \item{lower.kappa}{if \code{lower} is \code{NULL} this parameter is used to constrained optimization for the concentration parameter.}
  \item{upper.kappa}{if \code{upper} is \code{NULL} this parameter is used to constrained optimization for the concentration parameter.}
  \item{alpha}{if not \code{NULL} overrides the value of \code{p}. See the next argument \code{p}. This is a different
    parameterization, \code{alpha=-1/2} provides Hellinger distance,
    \code{alpha=-1} provides Kullback-Leibler distance and \code{alpha=-2}
    provides Neyman's Chi-Square distance.}
  \item{p}{\code{p=2} provides Hellinger distance, \code{p=-1}
	provides Kullback-Leibler distance and \code{p=Inf} provides Neyman's
	Chi-Square distance. It is ignored if \code{alpha} is not \code{NULL}.}
  \item{control.circular}{the attribute of the resulting object (\code{mu})}
  \item{digits}{integer indicating the precision to be used.}
  \item{\dots}{further parameters in \code{print.mde.vonmises}.}
}
\details{
The distance from an estimated density (by the non parametric kernel density estimator) and the model is evaluated by simple rectangular approximation. \code{optim} is used to performs minimization.
}

\value{
  Returns a list with the following components:

  \item{call}{the match.call().}
  \item{mu}{the estimate of the mean direction.}
  \item{kappa}{the estimate of the concentration parameter.}
  \item{dist}{the distance between the estimated density and the model.}
  \item{data}{the original supplied data converted in radians, clockwise and zero at 0.}
  \item{x}{the 'n' coordinates of the points where the density is estimated.}
  \item{y}{the estimated density values.}
  \item{k}{the density at the model.}
 }

 \references{
  C. Agostinelli. Robust estimation for circular data.
  Computational Statistics & Data Analysis, 51(12):5867-5875, 2007.
}

\author{Claudio Agostinelli}

\seealso{\code{\link[circular]{circular}}, \code{\link[circular]{mle.vonmises}} and \code{\link[wle]{wle.vonmises}}.}

\examples{

set.seed(1234)
x <- c(rvonmises(n=200, mu=circular(0), kappa=10), rvonmises(n=20, mu=circular(pi/2), kappa=20))
res <- mde.vonmises(x, bw=500, mu=circular(0), kappa=10)
res
plot(circular(0), type='n', xlim=c(-1, 1.75), shrink=1.2)
lines(circular(res$x), res$y)
lines(circular(res$x), res$k, col=2)
legend(1,1.5, legend=c('estimated density', 'MDE'), lty=c(1, 1), col=c(1, 2))

}

\keyword{robust}
