\name{plot.wle.cp}
\title{Plot the Weighted Mallows Cp}
\usage{
plot.wle.cp(object, base.line=0, num.max=20, 
           plot.it=TRUE, log.scale=FALSE, 
           xlab="Number of Predictors", ylab=NULL)
}

\alias{plot.wle.cp}

\arguments{
 \item{object}{an object of class \code{wle.cp}.}
 \item{base.line}{the intercept of the line to split the submodels in acceptable (good) and not-acceptable (bad), (the slope is always one).}
 \item{num.max}{maximum number of submodels plotted.}
 \item{plot.it}{if TRUE the graph is plotted.}
 \item{log.scale}{if TRUE the y-axis as log10 scale.}
 \item{xlab}{a title for the x axis.}
 \item{ylab}{a title for the y axis.}
}

\description{Plot the weighted Mallows Cp based on weighted likelihood.}

\value{
\item{num.good}{number of submodels below the \code{base.line}}
\item{num.bad}{number of submodels above the \code{base.line}}
\item{wcp.good}{list of the submodels below the \code{base.line} with their WCp.}
\item{wcp.bad}{list of the submodels above the \code{base.line} with their WCp.}
}


\references{

Agostinelli, C., (1999). Robust model selection in regression via weighted likelihood methodology, \emph{ Working Paper n. 1999.4}, Department of Statistics, Universiy of Padova.

Agostinelli, C., (1999). Robust model selection in regression via weighted likelihood methodology, submitted to \emph{ Statistics \& Probability Letters}, revised december 1999.

Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova.

Agostinelli, C., (1998). Verosimiglianza pesata nel modello di regressione lineare,  \emph{ XXXIX Riunione scientifica della Societ\`a Italiana di Statistica}, Sorrento 1998.

}

\author{Claudio Agostinelli}

%\note{}

\seealso{\code{\link{wle.cp}} a function to calculate the Weighted Mallows Cp, \code{\link{wle.lm}} a function for estimating linear models with normal distribution error and normal kernel.}

\examples{
library(wle)

x.data_c(runif(60,20,80),runif(5,73,78))
e.data_rnorm(65,0,0.6)
y.data_8*log(x.data+1)+e.data
y.data[61:65]_y.data[61:65]-4
z.data_c(rep(0,60),rep(1,5))

plot(x.data,y.data,xlab="X",ylab="Y")

xx.data <- cbind(x.data,x.data^2,x.data^3,log(x.data+1))

result <- wle.cp(y.data~xx.data)

plot(result,num.max=15)

}

\keyword{regression}
\keyword{robust}


