\name{mle.cp}
\title{Classical Mallows Cp}
\usage{
mle.cp(ydata, xdata, inter = 1, var.full = 0)
}
%- maybe also `usage' for other functions documented here.
\alias{mle.cp}
%- Also NEED an `\alias' for EACH other function documented here.
\arguments{
 \item{ydata}{A vector containing the observations from the response variable.}
 \item{xdata}{A matrix containing the observations from the explanatory variables.}
 \item{inter}{If 1, intercept is include in the model.}
 \item{var.full}{The value of variance to be used in the denominator of the Mallows Cp, if 0 the variance estimated from the full model is used.}
}
\description{The Mallows Cp is evaluated for each submodel.}
\value{
\item{cp}{Mallows Cp for each submodels}
\item{coefficients}{The parameters estimator, one row vector for eac submodel.}
\item{scale}{An estimation of the error scale, one value for each submodel.}
\item{residuals}{The residuals from the estimated model, one column vector for each submodel.}
\item{info}{Not well working yet, if 0 no error occurred.}
}

\author{Claudio Agostinelli}

%\seealso{}

\examples{
library(wle)

data(hald)

cor(hald)

result_mle.cp(y.hald,x.hald)

plot.mle.cp(result,num.max=6)
}

\keyword{Mallows Cp}
