\name{wle.stepwise}
\title{Weighted Stepwise, Backward and Forward selection methods}
\usage{
wle.stepwise(ydata,xdata,boot=100,group,inter=1,num.sol=1,raf=1,smooth=0.0320018,tol=10^(-6),equal=10^(-3),max.iter=500,min.weight=0.5,type="Forward",f.in=0.0,f.out=0.0,method="WLE")
}

%- maybe also `usage' for other functions documented here.
\alias{wle.stepwise}
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
 \item{min.weight}{See details.}
 \item{type}{
\code{type="Stepwise"}: The weighted stepwise methods is used,

\code{type="Forward"}: The weighted forward methods is used,

\code{type="Backward"}: The weighted backward method is used.}
 \item{f.in}{The in value}
 \item{f.out}{The out value}
 \item{method}{
\code{method="WLS"}: The submodel parameters are estimated by weighted least square with weights from the weighted likelihood estimator on the full model.

\code{method="WLE"}:The submodel parameters are estimated by weighted likelihood estimators.
}

}

\description{This function performs Weighted Stepwise, Forward and Backward model selection.}
\details{\code{min.weight}: the weighted likelihood equation could have more than one solution. These roots appear for particular situation depending on contamination level and type. The presence of multiple roots in the full model can create some problem in the set of weights we should use. Actualy, the selection of the root is done by the minimum scale error provided. Since this choice is not always the one would choose, we introduce the \code{min.weight} parameter in order to choose only between roots that do not down weight everything. This is not still the optimal solution, and perhaps, in the new release, this part will be change.}

\value{
\item{wstep}{the selected model}
\item{coefficients}{The parameters estimator, one row vector for each root found in the full model.}
\item{scale}{An estimation of the error scale, one value for each root found in the full model.}
\item{residuals}{The unweighted residuals from the estimated model, one column vector for each root found in the full model.}
\item{tot.weights}{The sum of the weights, one value for each root found in the full model.}
\item{weights}{The weights associated to each observation, one column vector for each root found in the full model.}
\item{tot.sol}{}
\item{freq}{The number of starting points converging to the roots.}
\item{index}{The root used as weighted likelihood estimator.}
\item{info}{Not well working yet, if 0 no error occurred.}
}
\references{

Agostinelli, C., (1999). Note on the Robust Stepwise Regression \emph{manuscript}.

Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova.

Agostinelli, C., (1998). Verosimiglianza pesata nel modello di regressione lineare,  \emph{ XXXIX Riunione scientifica della Societ\`a Italiana di Statistica}, Sorrento 1998.

}

\author{Claudio Agostinelli}

%\note{}

\seealso{\link{wle.smooth} an algorithm to choose the smoothing parameter for normal distribution and normal kernel, \link{wle.lm} a function for estimating linear models with normal distribution error and normal kernel.}

\examples{

library(wle)

# You can find this dataset in:
# Agostinelli, C., (1999). Robust model selection in regression 
# via weighted likelihood methodology, submitted to Statistics & 
# Probability Letters.

data(selection)

res_wle.stepwise(ydata,xdata,,boot=100,group=6,num.sol=3,min.weight=0.8,type="Stepwise",method="WLS",f.in=4,f.out=4)

res$wstep

plot(xdata[,1],ydata)

}

\keyword{Weighted Likelihood}
