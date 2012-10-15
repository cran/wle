\name{extractRoot}
\Rdversion{1.1}
\alias{extractRoot}
\alias{extractRoot.wle.glm}
\title{
Extract a Root from a result of a wle function
}
\description{
This function extract the information regarding one solution of the Weighted Likelihood Estimating Equation.

}
\usage{
\method{extractRoot}{wle.glm}(object, root=1, \dots)
}

\arguments{
  \item{object}{an object of class \code{"wle.glm"}, usually, a result of a
    call to \code{\link{wle.glm}}.}
  \item{root}{an integer number to specify which root should be extract.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
  \code{extract.wle.glm} returns an object of class \code{"extract.wle.glm.root"}, a (variable length) list 
  containing at least the following components:

    \item{coefficients}{a named vector of coefficients}
    \item{residuals}{the \emph{working} residuals, that is the residuals
    in the final iteration of the IWLS fit.  Since cases with zero
    weights are omitted, their working residuals are \code{NA}.}
    \item{fitted.values}{the fitted mean values, obtained by transforming
    the linear predictors by the inverse of the link function.}
    \item{rank}{the numeric rank of the fitted linear model.}
    \item{family}{the \code{\link{family}} object used.}
    \item{linear.predictors}{the linear fit on link scale.}
    \item{deviance}{up to a constant, minus twice the maximized
    log-likelihood.  Where sensible, the constant is chosen so that a
    saturated model has deviance zero.}
    \item{aic}{Akaike's \emph{An Information Criterion}, minus twice the
    maximized log-likelihood plus twice the number of coefficients (so
    assuming that the dispersion is known).}
    \item{null.deviance}{The deviance for the null model, comparable with
    \code{deviance}. The null model will include the offset, and an
    intercept if there is one in the model.  Note that this will be
    incorrect if the link function depends on the data other than
    through the fitted mean: specify a zero offset to force a correct
    calculation.}
    \item{iter}{the number of iterations of IWLS used.}
    \item{weights}{the \emph{working} weights, that is the weights
    in the final iteration of the IWLS fit.}
    \item{prior.weights}{the weights initially supplied, a vector of
    \code{1}s if none were.}
    \item{df.residual}{the residual degrees of freedom.}
    \item{df.null}{the residual degrees of freedom for the null model.}
    \item{y}{if requested (the default) the \code{y} vector
    used. (It is a vector even for a binomial model.)}
    \item{x}{if requested, the model matrix.}
    \item{model}{if requested (the default), the model frame.}
    \item{converged}{logical. Was the IWLS algorithm judged to have converged?}
    \item{boundary}{logical. Is the fitted value on the boundary of the
    attainable values?}
    \item{wle.weights}{final (robust) weights based on the WLE approach.}
    \item{wle.asymptotic}{logicals. If \code{TRUE} asymptotic weight based on Anscombe residual is used for the corresponding observation.}

    In addition, non-empty fits will have components \code{qr}, 
    \code{R}, \code{qraux}, \code{pivot}
    and \code{effects} relating to the final weighted linear fit.

  \item{family}{the \code{\link{family}} object used.}
  \item{call}{the matched call.}
  \item{formula}{the formula supplied.}
  \item{terms}{the \code{\link{terms}} object used.}
  \item{data}{the \code{data argument}.}
  \item{offset}{the offset vector used.}
  \item{control}{the value of the \code{control} argument used.}
  \item{method}{the name of the fitter function used, currently always
    \code{"wle.glm.fit"}.}
  \item{contrasts}{(where relevant) the contrasts used.}
  \item{xlevels}{(where relevant) a record of the levels of the factors
    used in fitting.}
  \item{tot.sol}{the number of solutions found.}
  \item{not.conv}{the number of starting points that does not converge after the \code{max.iter} (defined using \code{wle.glm.control}) iterations are reached.}
  \item{na.action}{(where relevant) information returned by
    \code{\link{model.frame}} on the special handling of \code{NA}s.}

  If a \code{\link{binomial}} \code{wle.glm} model was specified by 
  giving a two-column response, the weights returned by 
  \code{prior.weights} are
  the total numbers of cases (factored by the supplied case weights) and
  the component \code{y} of the result is the proportion of successes.
}

\author{
Claudio Agostinelli and Fatemah Al-quallaf
}

\seealso{
  \code{\link{anova.wle.glm.root}}
}

\examples{
## --- Continuing the Example from  '?wle.glm':
\dontshow{utils::example("wle.glm", echo = FALSE)}
anova(extractRoot(wle.glm.D93))
}

\keyword{utilities}

