\name{plot.mle.cp}
\title{Plot the Mallows Cp}
\usage{
plot.mle.cp(object, base.line=0, num.max=20, plot.it=TRUE, log.scale=FALSE)
}
%- maybe also `usage' for other functions documented here.
\alias{plot.mle.cp}
%- Also NEED an `\alias' for EACH other function documented here.
\arguments{
 \item{object}{An object from the function \code{mle.cp}.}
 \item{base.line}{The intercept of the line to split the submodels in acceptable (good) and not-acceptable (bad), (the slope is always 1).}
 \item{num.max}{Maximum number of submodels plotted.}
 \item{plot.it}{If TRUE the graph is plotted.}
 \item{log.scale}{If TRUE the y-axis as log10 scale.}
}
\description{Plot the Mallows Cp for each submodel.}
%\details{}

\value{
\item{num.good}{Number of submodels below the \code{base.line}}
\item{num.bad}{Number of submodels above the \code{base.line}}
\item{cp.good}{List of the submodels below the \code{base.line} with their Cp.}
\item{cp.bad}{List of the submodels above the \code{base.line} with their Cp.}
}

%\references{}

\author{Claudio Agostinelli}

%\note{}

%\seealso{}

\examples{
library(wle)

data(hald)

result_mle.cp(y.hald,x.hald)

plot.mle.cp(result,num.max=7)

}

\keyword{Mallows Cp}
