\name{wle.lm.summaries}
\title{Accessing Linear Model Fits for wle.lm}
\usage{
coef.wle.lm(object, ...)
formula.wle.lm(x, ...)
fitted.wle.lm(object, ...)
model.frame.wle.lm(formula, data, na.action, ...)
summary.wle.lm(object, root="ALL", ...)
summary.wle.lm.root(object, root=1, ...)

print.wle.lm(x, digits = max(3, getOption("digits") - 3), ...)

print.summary.wle.lm(x, digits = max(3, getOption("digits") - 3), 
           signif.stars= getOption("show.signif.stars"),  ...)

print.summary.wle.lm.root(x, digits = max(3, getOption("digits") - 3), 
           signif.stars= getOption("show.signif.stars"),  ...)
}

\alias{coef.wle.lm}
\alias{formula.wle.lm}
\alias{fitted.wle.lm}
\alias{model.frame.wle.lm}
\alias{summary.wle.lm}
\alias{summary.wle.lm.root}
\alias{weights.wle.lm}
\alias{print.wle.lm}
\alias{print.summary.wle.lm}
\alias{print.summary.wle.lm.root}

\arguments{
 \item{object}{an object of class \code{wle.lm}.}  
 \item{x}{an object of class \code{wle.lm} or \code{summary.wle.lm}.}
 \item{formula}{a model formula}
 \item{data}{\code{data.frame}, list, \code{environment} or object coercible to \code{data.frame} containing the variables in \code{formula}.}
 \item{na.action}{how \code{NA}s are treated. The default is first, any \code{na.action} attribute of \code{data}, second a \code{na.action} setting of \code{\link{options}}, and third \code{\link{na.fail}} if that is unset. The ``factory-fresh'' default is \code{\link{na.omit}}.}
 \item{root}{the root to be printed, in summary.wle.lm it could be  "ALL", all the roots are printed, or a vector of integers.}
 \item{digits}{number of digits to be used for most numbers.}
 \item{signif.stars}{logical; if \code{TRUE}, P-values are additionally encoded
   visually as ``significance stars'' in order to help scanning of long
   coefficient tables. It defaults to the \code{show.signif.stars} slot of \code{\link{options}}.}
 \item{...}{additional arguments.}
}
\description{All these functions are \code{\link{methods}} for class \code{wle.lm} or \code{summary.wle.lm}.}

\details{\code{print.summary.wle.lm} and \code{print.summary.wle.lm.root} tries formatting for each root the coefficients, standard errors, etc. and additionally gives ``significance stars'' if \code{signif.stars} is \code{TRUE}.
}

The generic accessor functions \code{coefficients}, \code{fitted.values}, \code{residuals} and \code{weights} can be used to extract various useful features of the value returned by \code{wle.lm}. 

\value{
  The function \code{summary.wle.lm} (the \code{summary.wle.lm.root} do the same for just one selected root) computes and returns, for each selected root, a list of summary statistics of the fitted linear model given in \code{object}, using the components (list elements) \code{"call"} and \code{"terms"} from its argument, plus 

 \item{residuals}{the \emph{weighted} residuals, the usual residuals rescaled by the square root of the weights given by \code{wle.lm}.}
 \item{coefficients}{a \eqn{p \times 4}{p x 4} matrix with columns for the estimated coefficient, its standard error, weighted-t-statistic and corresponding (two-sided) p-value.}
  \item{sigma}{the square root of the estimated variance of the random error.}
  \item{df}{degrees of freedom, a 3-vector \eqn{(p, \sum{weights} - p, p*)}.}
  \item{fstatistic}{a 3-vector with the value of the weighted-F-statistic with its numerator and denominator degrees of freedom.}
  \item{r.squared}{\eqn{R^2}, the ``fraction of variance explained by the model''.}
  \item{adj.r.squared}{the above \eqn{R^2} statistic
    \emph{``adjusted''}, penalizing for higher \eqn{p}.}
  \item{root}{the label of the root reported.}

}


\author{Claudio Agostinelli}

\seealso{\code{\link{wle.lm}} a function for estimating linear models with normal distribution error and normal kernel, \code{\link{plot.wle.lm}} for plot method.}

\examples{
library(wle)
# You can find this data set in:
# Hawkins, D.M., Bradu, D., and Kass, G.V. (1984). 
# Location of several outliers in multiple regression data using
# elemental sets. Technometrics, 26, 197-208.
#
data(artificial)

result <- wle.lm(y.artificial~x.artificial,boot=40,group=6,num.sol=3)

#summary only for the first root
summary(result,root=1)
#summary for all the roots
summary(result,root="ALL")
}


\keyword{robust}
\keyword{regression}

