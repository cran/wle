\name{plot.wle.cp}
\title{Plot the Weighted Mallows Cp}
\usage{
plot.wle.cp(object, base.line=0, num.max=20, plot.it=TRUE, log.scale=FALSE)
}
%- maybe also `usage' for other functions documented here.
\alias{plot.wle.cp}
%- Also NEED an `\alias' for EACH other function documented here.
\arguments{
 \item{object}{An object from the function \code{wle.cp}.}
 \item{base.line}{The intercept of the line to split the submodels in acceptable (good) and not-acceptable (bad), (the slope is always 1).}
 \item{num.max}{Maximum number of submodels plotted.}
 \item{plot.it}{If TRUE the graph is plotted.}
 \item{log.scale}{If TRUE the y-axis as log10 scale.}
}
\description{Plot the WCp for each submodel.}
%\details{}

\value{
\item{num.good}{Number of submodels below the \code{base.line}}
\item{num.bad}{Number of submodels above the \code{base.line}}
\item{wcp.good}{List of the submodels below the \code{base.line} with their WCp.}
\item{wcp.bad}{List of the submodels above the \code{base.line} with their WCp.}
}

\references{

Agostinelli, C., (1999). Robust model selection in regression via weighted likelihood methodology, \emph{ Working Paper n. 1999.4}, Department of Statistics, Universiy of Padova.

Agostinelli, C., (1999). Robust model selection in regression via weighted likelihood methodology, submitted to \emph{ Statistics \& Probability Letters}.

Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova.

Agostinelli, C., (1998). Verosimiglianza pesata nel modello di regressione lineare,  \emph{ XXXIX Riunione scientifica della Societ\`a Italiana di Statistica}, Sorrento 1998.

}

\author{Claudio Agostinelli}

%\note{}

\seealso{\link{wle.cp} a function to calculate the Weighted Mallows Cp, \link{wle.lm} a function for estimating linear models with normal distribution error and normal kernel.}

\examples{
library(wle)

x.data_c(runif(60,20,80),runif(5,73,78))
e.data_rnorm(65,0,0.6)
y.data_8*log(x.data+1)+e.data
y.data[61:65]_y.data[61:65]-4
z.data_c(rep(0,60),rep(1,5))

plot(x.data,y.data,xlab="X",ylab="Y")

xx.data_cbind(x.data,x.data^2,x.data^3,log(x.data+1))

result_wle.cp(y.data,xx.data,boot=20,group=10,num.sol=3)

plot.wle.cp(result,num.max=15)

}

\keyword{Weighted Likelihood}
