\name{mle.cv}
\title{Cross Validation Selection Method}
\usage{
mle.aic(ydata, xdata, monte.carlo = 500, split, inter = 1)
}
%- maybe also `usage' for other functions documented here.
\alias{mle.cv}
%- Also NEED an `\alias' for EACH other function documented here.
\arguments{
 \item{ydata}{A vector containing the observations from the response variable.}
 \item{xdata}{A matrix containing the observations from the explanatory variables.}
 \item{monte.carlo}{The number of Monte Carlo replication we use to estimate the average prediction error.}
 \item{split}{The size of the costruction sample. When the suggested value is outside the possible range, the split size is let equal to \eqn{max(round(size^{(3/4)}),nvar+inter+2)}.}
 \item{inter}{If 1, intercept is include in the model.}
}
\description{The Cross Validation selection method is evaluated for each submodel.}
\value{
\item{cv}{The estimated error for each submodels}
\item{info}{Not well working yet, if 0 no error occurred.}
}

\author{Claudio Agostinelli}

%\seealso{}

\examples{
library(wle)

data(hald)

cor(hald)

result_mle.cv(y.hald,x.hald,split=0)

result$cv[order(result$cv[,6]),][1:10,]
}

\keyword{Cross Validation}
