\name{wle.cp}
\title{Weighted Mallows Cp}

\usage{
wle.cp(formula, data=list(), model=TRUE, x=FALSE, 
       y=FALSE, boot=30, group, var.full=0, num.sol=1, 
       raf="HD", smooth=0.031, tol=10^(-6), 
       equal=10^(-3), max.iter=500, min.weight=0.5, 
       method="full", alpha=2, contrasts=NULL, verbose=FALSE)
}
\alias{wle.cp}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.
    The details of model specification are given below.}

  \item{data}{an optional data frame containing the variables
    in the model.  By default the variables are taken from
    the environment which \code{wle.cp} is called from.}

  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding components of the fit (the model frame, the model matrix, the
    response.)}

  \item{boot}{the number of starting points based on boostrap subsamples to use in the search of the roots.}
  
  \item{group}{the dimension of the bootstap subsamples. The default value is \eqn{max(round(size/4),var)} where \eqn{size} is the number of observations and \eqn{var} is the number of variables.}

  \item{var.full}{the value of variance to be used in the denominator of the WCP, if 0 the variance estimated from the full model is used.}

  \item{num.sol}{maximum number of roots to be searched.}

  \item{raf}{type of Residual adjustment function to be use:

\code{raf="HD"}: Hellinger Distance RAF,

\code{raf="NED"}: Negative Exponential Disparity RAF,

\code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}

  \item{smooth}{the value of the smoothing parameter.}

  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}c

  \item{equal}{the absolute value for which two roots are considered the same. (This parameter must be greater than \code{tol}).}

  \item{max.iter}{maximum number of iterations.}

  \item{min.weight}{see details.}

  \item{method}{see details.}

  \item{alpha}{penalty value.}

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}

  \item{verbose}{if TRUE warnings are printed.}
}

\description{The Weighted Mallows Cp is evaluated for each submodel.}

\details{
Models for \code{wle.cp} are specified symbolically.  A typical model has the form \code{response ~ terms} where \code{response} is the (numeric) response vector and \code{terms} is a series of terms which specifies a linear predictor for \code{response}.  A terms specification of the form \code{first+second} indicates all the terms in \code{first} together with all the terms in \code{second} with duplicates removed. A specification of the form \code{first:second} indicates the the set of terms obtained by taking the interactions of all terms in \code{first} with all terms in \code{second}. The specification \code{first*second} indicates the \emph{cross} of \code{first} and \code{second}.  This is the same as \code{first+second+first:second}.

\code{min.weight}: the weighted likelihood equation could have more than one solution. These roots appear for particular situation depending on contamination level and type. The presence of multiple roots in the full model can create some problem in the set of weights we should use. Actually, the selection of the root is done by the minimum scale error provided. Since this choice is not always the one would choose, we introduce the \code{min.weight} parameter in order to choose only between roots that do not down weight everything. This is not still the optimal solution, and perhaps, in the new release, this part will be change. 

\code{method}: this parameter, when set to "reduced",  allows to use weights based on the reduced model. This is strongly discourage since the robust and asymptotic property of this kind of weighted Cp are not as good as the one based on \code{method="full"}.}

\value{

\code{wle.cp} returns an object of \code{\link{class}} \code{"wle.cp"}.

  The function \code{summary} is used to obtain and print a summary of the results.
  The generic accessor functions \code{coefficients} and \code{residuals} extract coefficients and residuals returned by \code{wle.cp}.
  The object returned by \code{wle.cp} are:

\item{wcp}{Weighted Mallows Cp for each submodels}
\item{coefficients}{the parameters estimator, one row vector for each root found and each submodel.}
\item{scale}{an estimation of the error scale, one value for each root found and each submodel.}
\item{residuals}{the unweighted residuals from the estimated model, one column vector for each root found and each submodel.}
\item{tot.weights}{the sum of the weights divide by the number of observations, one value for each root found and each submodel.}
\item{weights}{the weights associated to each observation, one column vector for each root found and each submodel.}
\item{freq}{the number of starting points converging to the roots.}
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

Agostinelli, C., (1999). Robust model selection in regression via weighted likelihood methodology, \emph{ Working Paper n. 1999.4}, Department of Statistics, Universiy of Padova.

Agostinelli, C., (1999). Robust model selection in regression via weighted likelihood methodology, submitted to \emph{ Statistics \& Probability Letters}, revised december 1999.

Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova.

Agostinelli, C., Markatou, M., (1998). A one-step robust estimator for regression based on the weighted likelihood reweighting scheme, \emph{ Statistics \& Probability Letters}, Vol. 37, n. 4, 341-350.

Agostinelli, C., (1998). Verosimiglianza pesata nel modello di regressione lineare,  \emph{ XXXIX Riunione scientifica della Societ\`a Italiana di Statistica}, Sorrento 1998.

}

\author{Claudio Agostinelli}

\seealso{\link{wle.smooth} an algorithm to choose the smoothing parameter for normal distribution and normal kernel, \link{wle.lm} a function for estimating linear models with normal distribution error and normal kernel.}

\examples{
library(wle)

x.data <- c(runif(60,20,80),runif(5,73,78))
e.data <- rnorm(65,0,0.6)
y.data <- 8*log(x.data+1)+e.data
y.data[61:65] <- y.data[61:65]-4
z.data <- c(rep(0,60),rep(1,5))

plot(x.data,y.data,xlab="X",ylab="Y")

xx.data <- cbind(x.data,x.data^2,x.data^3,log(x.data+1))
colnames(xx.data) <- c("X","X^2","X^3","log(X+1)")

result <- wle.cp(y.data~xx.data,boot=10,group=10,num.sol=2)

summary(result)

plot(result,num.max=15)

result <- wle.cp(y.data~xx.data+z.data,boot=10,group=10,num.sol=2)

summary(result)

plot(result,num.max=15)
}

\keyword{regression}
\keyword{robust}






