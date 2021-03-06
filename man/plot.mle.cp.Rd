\name{plot.mle.cp}
\title{Plot the Mallows Cp}
\usage{
\method{plot}{mle.cp}(x, base.line=0, num.max=20, 
           plot.it=TRUE, log.scale=FALSE, 
           xlab="Number of Predictors", ylab=NULL,
           verbose=FALSE, ...)
}

\alias{plot.mle.cp}

\arguments{
 \item{x}{an object of class \code{mle.cp}.}
 \item{base.line}{the intercept of the line to split the submodels in acceptable (good) and not-acceptable (bad), (the slope is always one).}
 \item{num.max}{maximum number of submodels plotted.}
 \item{plot.it}{if TRUE the graph is plotted.}
 \item{log.scale}{if TRUE the y-axis as log10 scale.}
 \item{xlab}{a title for the x axis.}
 \item{ylab}{a title for the y axis.}
 \item{verbose}{if \code{TRUE} warnings are printed.}
 \item{...}{graphical parameters can be given as arguments.}
}

\description{Plot the Mallows Cp.}

\value{
\item{num.good}{number of submodels below the \code{base.line}}
\item{num.bad}{number of submodels above the \code{base.line}}
\item{cp.good}{list of the submodels below the \code{base.line} with their Cp.}
\item{cp.bad}{list of the submodels above the \code{base.line} with their Cp.}
}

%\references{}

\author{Claudio Agostinelli}

%\note{}

\seealso{\code{\link{mle.cp}} a function to calculate the Mallows Cp.}

\examples{
library(wle)

data(hald)

result <- mle.cp(y.hald~x.hald)

plot(result,num.max=7)

}

\keyword{regression}


