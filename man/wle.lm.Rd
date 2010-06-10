\name{wle.lm}
\title{Fitting Linear Models using Weighted Likelihood}
\usage{
wle.lm(formula, data=list(), model=TRUE, x=FALSE, 
       y=FALSE, boot=30, group, num.sol=1, raf="HD", 
       smooth=0.031, tol=10^(-6), equal=10^(-3), 
       max.iter=500, contrasts=NULL, verbose=FALSE)
}

\alias{wle.lm}

\arguments{
  \item{formula}{a symbolic description of the model to be fit.
    The details of model specification are given below.}

  \item{data}{an optional data frame containing the variables
    in the model.  By default the variables are taken from
    the environment which \code{wle.lm} is called from.}

  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding components of the fit (the model frame, the model matrix, the
    response.)}

  \item{boot}{the number of starting points based on boostrap subsamples to use in the search of the roots.}

  \item{group}{the dimension of the bootstap subsamples. The default value is \eqn{max(round(size/4),var)} where \eqn{size} is the number of observations and \eqn{var} is the number of variables.}

  \item{num.sol}{maximum number of roots to be searched.}

  \item{raf}{type of Residual adjustment function to be used:

\code{raf="HD"}: Hellinger Distance RAF,

\code{raf="NED"}: Negative Exponential Disparity RAF,

\code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}

  \item{smooth}{the value of the smoothing parameter.}

  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}

  \item{equal}{the absolute value for which two roots are considered the same. (This parameter must be greater than \code{tol}).}

  \item{max.iter}{maximum number of iterations.}

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}

  \item{verbose}{if \code{TRUE} warnings are printed.}
}
\description{
\code{wle.lm} is used to fit linear models via Weighted Likelihood, when the errors are iid from a normal distribution with null mean and unknown variance. The carriers are considered fixed. Note that this estimator is robust against the presence of bad leverage points too.
}

\details{

Models for \code{wle.lm} are specified symbolically.  A typical model has the form \code{response ~ terms} where \code{response} is the (numeric) response vector and \code{terms} is a series of terms which specifies a linear predictor for \code{response}.  A terms specification of the form \code{first+second} indicates all the terms in \code{first} together with all the terms in \code{second} with duplicates removed. A specification of the form \code{first:second} indicates the the set of terms obtained by taking the interactions of all terms in \code{first} with all terms in \code{second}. The specification \code{first*second} indicates the \emph{cross} of \code{first} and \code{second}.  This is the same as \code{first+second+first:second}.
}

\value{

\code{wle.lm} returns an object of \code{\link{class}} \code{"wle.lm"}.

  The function \code{summary} is used to obtain and print a summary of the results.
  The generic accessor functions \code{coefficients}, \code{residuals} and \code{fitted.values} extract coefficients, residuals and fitted values returned by \code{wle.lm}.

  The object returned by \code{wle.lm} are:

\item{coefficients}{the parameters estimator, one row vector for each root found.}
\item{standard.error}{an estimation of the standard error of the parameters estimator, one row vector for each root found.}
\item{scale}{an estimation of the error scale, one value for each root found.}
\item{residuals}{the unweighted residuals from the estimated model, one column vector for each root found.}
\item{fitted.values}{the fitted values from the estimated model, one column vector for each root found.}
\item{tot.weights}{the sum of the weights divide by the number of observations, one value for each root found.}
\item{weights}{the weights associated to each observation, one column vector for each root found.}
\item{f.density}{the non-parametric density estimation.}
\item{m.density}{the smoothed model.}
\item{delta}{the Pearson residuals.}
\item{freq}{the number of starting points converging to the roots.}
\item{tot.sol}{the number of solutions found.}
\item{not.conv}{the number of starting points that does not converge after the \code{max.iter} iterations are reached.}
\item{call}{the match.call().}
\item{contrasts}{}
\item{xlevels}{}
\item{terms}{the model frame.}
\item{model}{if \code{model=TRUE} a matrix with first column the dependent variable and the remain column the explanatory variables for the full model.}
\item{x}{if \code{x=TRUE} a matrix with the explanatory variables for the full model.}
\item{y}{if \code{y=TRUE} a vector with the dependent variable.}
\item{info}{not well working yet, if 0 no error occurred.}
}

\references{
Agostinelli, C., (1998) Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{Ph.D Thesis}, Department of Statistics, University of Padova.

Agostinelli, C., Markatou, M., (1998) A one-step robust estimator for regression based on the weighted likelihood reweighting scheme, \emph{ Statistics \& Probability Letters}, Vol. 37, n. 4, 341-350.

Agostinelli, C., (1998) Verosimiglianza pesata nel modello di regressione lineare,  \emph{ XXXIX Riunione scientifica della Societ\`a Italiana di Statistica}, Sorrento 1998.

}

\author{Claudio Agostinelli}

%\note{}

\seealso{\link{wle.smooth} an algorithm to choose the smoothing parameter for normal distribution and normal kernel.}

\examples{
library(wle)
# You can find this data set in:
# Hawkins, D.M., Bradu, D., and Kass, G.V. (1984). 
# Location of several outliers in multiple regression data using
# elemental sets. Technometrics, 26, 197-208.
#
data(artificial)

result <- wle.lm(y.artificial~x.artificial,boot=40,num.sol=3)

summary(result)

plot(result)
}

\keyword{robust}
\keyword{regression}


