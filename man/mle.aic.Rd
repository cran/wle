\name{mle.aic}
\title{Classical Akaike Information Criterion}
\usage{
mle.aic(ydata, xdata, inter = 1, var.full = 0, alpha = 2)
}
%- maybe also `usage' for other functions documented here.
\alias{mle.aic}
%- Also NEED an `\alias' for EACH other function documented here.
\arguments{
 \item{ydata}{A vector containing the observations from the response variable.}
 \item{xdata}{A matrix containing the observations from the explanatory variables.}
 \item{inter}{If 1, intercept is include in the model.}
 \item{var.full}{The value of variance to be used, if 0 the variance estimated from the full model is used.}
 \item{alpha}{The penalized constant.}
}
\description{The Akaike Information Criterion is evaluated for each submodel.}
\value{
\item{aic}{The AIC for each submodels}
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

result_mle.aic(y.hald,x.hald)

result$aic[order(result$aic[,6]),][1:10,]
}

\keyword{Akaike Information Criterion, AIC}
