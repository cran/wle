\name{anova.wle.glm.root}
\alias{anova.wle.glm.root}
\title{Robust Analysis of Deviance for Generalized Linear Model Fits}
\usage{
\method{anova}{wle.glm.root}(object, \dots, dispersion = NULL, test = NULL)
}
\alias{anova.wleglmlist}
\description{
  Compute a robust analysis of deviance table for one or more generalized
  linear model fits.
}
\arguments{
  \item{object, \dots}{objects of class \code{wle.glm.root}, typically
    the result of a call to \code{\link{extractRoot.wle.glm}}, or a list of
    \code{objects} each of which a result of a call to \code{"extractRoot.wle.glm"} method.}
  \item{dispersion}{the dispersion parameter for the fitting family.
    By default it is obtained from the object(s).}
  \item{test}{a character string, (partially) matching one of \code{"Chisq"},
    \code{"F"} or \code{"Cp"}. See \code{\link{stat.anova}}.}
}
\details{
  Specifying a single object gives a sequential analysis of deviance
  table for that fit.  That is, the reductions in the residual deviance
  as each term of the formula is added in turn are given in as
  the rows of a table, plus the residual deviances themselves.

  If more than one object is specified, the table has a row for the
  residual degrees of freedom and deviance for each model.
  For all but the first model, the change in degrees of freedom and
  deviance is also given. (This only makes statistical sense if the
  models are nested.)  It is conventional to list the models from
  smallest to largest, but this is up to the user.

  The table will optionally contain test statistics (and P values)
  comparing the reduction in deviance for the row to the residuals.
  For models with known dispersion (e.g., binomial and Poisson fits)
  the robust chi-squared test is most appropriate, and for those with
  dispersion estimated by moments (e.g., \code{gaussian},
  \code{quasibinomial} and \code{quasipoisson} fits) the Robust F test is
  most appropriate.  Robust Mallows' \eqn{C_p}{Cp} statistic is the residual
  weighted deviance plus twice the estimate of \eqn{\sigma^2}{sigma^2} times
  the residual (weighted) degrees of freedom, which is closely related to 
  Robust AIC (and a multiple of it if the dispersion is known).

  The dispersion estimate will be taken from the largest model, using
  the value returned by \code{\link{summary.wle.glm}}. As this will in most
  cases use a Chisquared-based estimate, the F tests are not based on
  the residual deviance in the analysis of deviance table shown.
}
\value{
  An object of class \code{"anova"} inheriting from class \code{"data.frame"}.
}
\section{Warning}{
  The comparison between two or more models by \code{anova.wle.glm.root} or
  \code{anova.wleglmlist} will only be valid if they
  are fitted to the same dataset. This may be a problem if there are
  missing values and \R's default of \code{na.action = na.omit} is used,
  and \code{anova.wleglmlist} will detect this with an error.

  Since in a model selection procedure and/or on an ANOVA table the weights of the WLE procedure must be that of the FULL model (and not that of the actual model) statistics on degrees of freedom, deviance and AIC are valid only if \code{object} is the FULL model.
}
\references{
  Agostinelli, C. and Markatou, M. (2001) Test of hypotheses based on 
  the Weighted Likelihood Methodology, \emph{Statistica Sinica},  
  vol. 11, n. 2, 499-514.

  Agostinelli, C. (2002) Robust model selection in regression via weighted 
  likelihood methodology Statistics and Probability Letters, 56, 289-300.

  Agostinelli, C. and Al-quallaf, F. (2009) Robust inference in 
  Generalized Linear Models. Manuscript in preparation.

  Hastie, T. J. and Pregibon, D. (1992)
  \emph{Generalized linear models.}
  Chapter 6 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}
\seealso{
  \code{\link{extractRoot.wle.glm}}, \code{\link{wle.glm}}, \code{\link{anova}}.
}
\examples{
## --- Continuing the Example from  '?wle.glm':
\dontshow{require(utils)
example("wle.glm", echo = FALSE)}
anova(extractRoot(wle.glm.D93))
anova(extractRoot(wle.glm.D93), test = "Cp")
anova(extractRoot(wle.glm.D93), test = "Chisq")
}

\keyword{robust}
\keyword{models}
\keyword{regression}

