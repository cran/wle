\name{wle.lm}
\title{Fitting Linear Models}
\usage{
wle.lm(ydata, xdata, boot=100, group, inter=1, num.sol=1, raf=1, smooth=0.0320018, tol=10^(-6), equal=10^(-3), max.iter=500)
}
%- maybe also `usage' for other functions documented here.
\alias{wle.lm}
%- Also NEED an `\alias' for EACH other function documented here.
\arguments{
 \item{ydata}{A vector containing the observations from the response variable.}
 \item{xdata}{A matrix containing the observations from the explanatory variables.}
 \item{boot}{The number of starting points based on boostrap subsamples to use in the search of the roots.}
 \item{group}{The dimension of the bootstap subsamples.}
 \item{inter}{If 1, intercept is include in the model.}
 \item{num.sol}{Maximum number of roots to be searched.}
 \item{raf}{Type of Residual adjustment function to be use:

\code{raf=1}: Hellinger Distance RAF,

\code{raf=2}: Negative Exponential Disparity RAF,

\code{raf=3}: Chi-Squared Disparity RAF.}
 \item{smooth}{The value of the smoothing parameter.}
 \item{tol}{The absolute accuracy to be used to achieve convergence of the algorithm.}
 \item{equal}{The absolute value for which two roots are considered the same. (This parameter must be greater than \code{tol}).}
 \item{max.iter}{Maximum number of iterations.}
}
\description{
\code{wle.lm} is used to fit linear models via Weighted Likelihood Methodology, when the errors are iid from a normal distribution with null mean and unknown variance. The carriers are considered fixed and no particular attention is devoted to the presence of leverage points. It will developed a version to deals with the presence of leverage points.
}
%\details{}
\value{
\item{coefficients}{The parameters estimator, one row vector for each root found.}

\item{standard.error}{An estimation of the standard error of the parameters estimator, one row vector for each root found.}

\item{scale}{An estimation of the error scale, one value for each root found.}

\item{residuals}{The unweighted residuals from the estimated model, one column vector for each root found.}

\item{tot.weights}{The sum of the weights, one value for each root found.}

\item{weights}{The weights associated to each observation, one column vector for each root found.}

\item{freq}{The number of starting points converging to the roots.}

\item{tot.sol}{The number of solutions found.}

\item{not.conv}{The number of starting points that does not converge after the \code{imax} iteration are reached.}

}

\references{
Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova.

Agostinelli, C., Markatou, M., (1998). A one-step robust estimator for regression based on the weighted likelihood reweighting scheme, \emph{ Statistics \& Probability Letters}, Vol. 37, n. 4, 341-350.

Agostinelli, C., (1998). Verosimiglianza pesata nel modello di regressione lineare,  \emph{ XXXIX Riunione scientifica della Societ\`a Italiana di Statistica}, Sorrento 1998, in print.

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
result_wle.lm(y.artificial,x.artificial,boot=40,group=6,num.sol=3)
print(result$coeff)
plot.wle.lm(result)
}

\keyword{Weighted Likelihood}
