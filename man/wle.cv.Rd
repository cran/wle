\name{wle.cv}
\title{Model Selection by Weighted Cross-Validation}
\usage{
wle.cv(ydata, xdata, boot=100, group, monte.carlo=500, split, inter=1, num.sol=1, raf=1, smooth=0.0320018, tol=10^(-6), equal=10^(-3), max.iter=500, min.weight=0.5)
}
\alias{wle.cv}
\arguments{
 \item{ydata}{A vector containing the observations from the response variable.}
 \item{xdata}{A matrix containing the observations from the explanatory variables.}
 \item{boot}{The number of starting points based on boostrap subsamples to use in the search of the roots.}
 \item{group}{The dimension of the bootstap subsamples.}
 \item{monte.carlo}{The number of Monte Carlo replication we use to estimate the average prediction error.}
 \item{split}{The size of the costruction sample. When the suggested value is outside the possible range, the split size is let equal to \eqn{max(round(size^{(3/4)}),nvar+inter+2)}.}
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
 \item{min.weight}{See details.}
}
\description{
The weighted Cross-Validation methods is used to choose the best linear model.
}
\details{
\code{min.weight}: the weighted likelihood equation could have more than one solution. These roots appear for particular situation depending on contamination level and type. The presence of multiple roots in the full model can create some problem in the set of weights we should use. Actualy, the selection of the root is done by the minimum scale error provided. Since this choice is not always the one would choose, we introduce the \code{min.weight} parameter in order to choose only between roots that do not down weight everything. This is not still the optimal solution, and perhaps, in the new release, this part will be change.}

\value{
\item{wcv}{Weighted Cross-Validation for each submodels}
\item{coefficients}{The parameters estimator, one row vector for each root found and each submodel.}
\item{scale}{An estimation of the error scale, one value for each root found and each submodel.}
\item{residuals}{The unweighted residuals from the estimated model, one column vector for each root found and each submodel.}
\item{tot.weights}{The sum of the weights, one value for each root found and each submodel.}
\item{weights}{The weights associated to each observation, one column vector for each root found and each submodel.}
\item{freq}{The number of starting points converging to the roots.}
\item{index}{Position of the root used for the weights.}
\item{info}{Not well working yet, if 0 no error occurred.}
}

\references{
Agostinelli, C., (1999). Robust model selection by Cross-Validation via weighted likelihood methodology, \emph{Working Paper n. 1999.37}, Department of Statistics, Universiy of Padova.

Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova.

Agostinelli, C., Markatou, M., (1998). A one-step robust estimator for regression based on the weighted likelihood reweighting scheme, \emph{ Statistics \& Probability Letters}, Vol. 37, n. 4, 341-350.

Agostinelli, C., (1998). Verosimiglianza pesata nel modello di regressione lineare,  \emph{ XXXIX Riunione scientifica della Societ\`a Italiana di Statistica}, Sorrento 1998.
}

\author{Claudio Agostinelli}

%\note{}

\seealso{\link{wle.smooth} an algorithm to choose the smoothing parameter for normal distribution and normal kernel, \link{wle.lm} a function for estimating linear models with normal distribution error and normal kernel.}

\examples{
library(wle)

x.data_c(runif(60,20,80),runif(5,73,78))
e.data_rnorm(65,0,0.6)
y.data_8*log(x.data+1)+e.data
y.data[61:65]_y.data[61:65]-4

plot(x.data,y.data,xlab="X",ylab="Y")

xx.data_cbind(x.data,x.data^2,x.data^3,log(x.data+1))

result_wle.cv(y.data,xx.data,boot=20,group=10,num.sol=3,split=0)

result$wcv[order(result$wcv[,6]),][1:10,]
}

\keyword{Weighted Likelihood}