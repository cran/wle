\name{wle.stepwise}
\title{Weighted Stepwise, Backward and Forward selection methods}
\usage{
wle.stepwise(formula, data=list(), model=TRUE, x=FALSE, 
             y=FALSE, boot=30, group, num.sol=1, raf="HD", 
             smooth=0.031, tol=10^(-6), equal=10^(-3), 
             max.iter=500, min.weight=0.5, type="Forward", 
             f.in=4.0, f.out=4.0, method="WLE", 
             contrasts=NULL, verbose=FALSE)
}
\alias{wle.stepwise}

\arguments{

  \item{formula}{a symbolic description of the model to be fit.
    The details of model specification are given below.}

  \item{data}{an optional data frame containing the variables
    in the model.  By default the variables are taken from
    the environment which \code{wle.stepwise} is called from.}

  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding components of the fit (the model frame, the model matrix, the
    response.)}

  \item{boot}{the number of starting points based on boostrap subsamples to use in the search of the roots.}

  \item{group}{the dimension of the bootstap subsamples. The default value is \eqn{max(round(size/4),var)} where \eqn{size} is the number of observations and \eqn{var} is the number of variables.}

  \item{num.sol}{maximum number of roots to be searched.}

  \item{raf}{type of Residual adjustment function to be use:

\code{raf="HD"}: Hellinger Distance RAF,

\code{raf="NED"}: Negative Exponential Disparity RAF,

\code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}
 
  \item{smooth}{the value of the smoothing parameter.}

  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}

  \item{equal}{the absolute value for which two roots are considered the same. (This parameter must be greater than \code{tol}).}

  \item{max.iter}{maximum number of iterations.}

  \item{min.weight}{see details.}

  \item{type}{
\code{type="Stepwise"}: the weighted stepwise methods is used,

\code{type="Forward"}: the weighted forward methods is used,

\code{type="Backward"}: the weighted backward method is used.}

  \item{f.in}{the in value}

  \item{f.out}{the out value}

  \item{method}{
\code{method="WLS"}: the submodel parameters are estimated by weighted least square with weights from the weighted likelihood estimator on the full model.

\code{method="WLE"}: the submodel parameters are estimated by weighted likelihood estimators.
}

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}

  \item{verbose}{if \code{TRUE} warnings are printed.}
}

\description{This function performs Weighted Stepwise, Forward and Backward model selection.}

\details{

Models for \code{wle.stepwise} are specified symbolically.  A typical model has the form \code{response ~ terms} where \code{response} is the (numeric) response vector and \code{terms} is a series of terms which specifies a linear predictor for \code{response}.  A terms specification of the form \code{first+second} indicates all the terms in \code{first} together with all the terms in \code{second} with duplicates removed. A specification of the form \code{first:second} indicates the the set of terms obtained by taking the interactions of all terms in \code{first} with all terms in \code{second}. The specification \code{first*second} indicates the \emph{cross} of \code{first} and \code{second}.  This is the same as \code{first+second+first:second}.

\code{min.weight}: the weighted likelihood equation could have more than one solution. These roots appear for particular situation depending on contamination level and type. The presence of multiple roots in the full model can create some problem in the set of weights we should use. Actually, the selection of the root is done by the minimum scale error provided. Since this choice is not always the one would choose, we introduce the \code{min.weight} parameter in order to choose only between roots that do not down weight everything. This is not still the optimal solution, and perhaps, in the new release, this part will be change. }

\value{


\code{wle.stepwise} returns an object of \code{\link{class}} \code{"wle.stepwise"}.

  The function \code{summary} is used to obtain and print a summary of the results.
  The generic accessor functions \code{coefficients} and \code{residuals} extract coefficients and residuals returned by \code{wle.stepwise}.

  The object returned by \code{wle.stepwise} are:

\item{wstep}{the iterations with the model selected.}
\item{coefficients}{the parameters estimator, one row vector for each root found in the full model.}
\item{scale}{an estimation of the error scale, one value for each root found in the full model.}
\item{residuals}{the unweighted residuals from the estimated model, one column vector for each root found in the full model.}
\item{tot.weights}{the sum of the weights divide by the number of observations, one value for each root found in the full model.}
\item{weights}{the weights associated to each observation, one column vector for each root found in the full model.}
\item{freq}{the number of starting points converging to the roots.}
\item{index}{position of the root used for the weights.}
\item{call}{the match.call().}
\item{contrasts}{}
\item{xlevels}{}
\item{terms}{the model frame.}
\item{model}{if \code{model=TRUE} a matrix with first column the dependent variable and the remain column the explanatory variables for the full model.}
\item{x}{if \code{x=TRUE} a matrix with the explanatory variables for the full model.}
\item{y}{if \code{y=TRUE} a vector with the dependent variable.}
\item{info}{not well working yet, if 0 no error occurred.}
  \item{type}{\code{"Stepwise"}: the weighted stepwise methods is used, \code{"Forward"}: the weighted forward methods is used, \code{"Backward"}: the weighted backward method is used.}

  \item{f.in}{the in value.}

  \item{f.out}{the out value.}

  \item{method}{if "WLS" the submodel parameters are estimated by weighted least square with weights from the weighted likelihood estimator on the full model else if "WLE" the submodel parameters are estimated by weighted likelihood estimators.
}

}

\references{

Agostinelli, C., (2000) Robust stepwise regression, Working Paper n. 2000.10 del Dipartimento di Scienze Statistiche, Universit\`a di Padova, Padova.

Agostinelli, C., (2001) Robust stepwise regression, \emph{Journal of Applied Statistics} in press.

Agostinelli, C., (1998) Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{Ph.D Thesis}, Department of Statistics, University of Padova.

Agostinelli, C., (1998) Verosimiglianza pesata nel modello di regressione lineare,  \emph{XXXIX Riunione scientifica della Societ\`a Italiana di Statistica}, Sorrento 1998.

}

\author{Claudio Agostinelli}

\seealso{\link{wle.smooth} an algorithm to choose the smoothing parameter for normal distribution and normal kernel, \link{wle.lm} a function for estimating linear models with normal distribution error and normal kernel.}

\examples{

library(wle)

# You can find this dataset in:
# Agostinelli, C., (1999). Robust model selection in regression 
# via weighted likelihood methodology, submitted to Statistics & 
# Probability Letters.

data(selection)

result <- wle.stepwise(ydata~xdata,,boot=100,group=6,
            num.sol=3,min.weight=0.8,type="Stepwise",
            method="WLS")

summary(result)
}

\keyword{robust}
\keyword{regression}





