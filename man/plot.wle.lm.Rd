\name{plot.wle.lm}
\title{Plots for the Linear Model}
\usage{
plot.wle.lm(object, level.weight=0.5, plot.it=TRUE)
}
\alias{plot.wle.lm}

\arguments{
 \item{object}{Object from the \code{wle.lm} function output.}
 \item{level.weight}{Value of the weight under which an observations is marked with different color.}
 \item{plot.it}{Not working yet.}
}
\description{The \code{plot.wle.lm} function plot a separate graph windows for each roots. In each windows four plots are printed: residuals vs fitted, normal qq plot of the residuals, weighted residuals vs weighted fitted, normal qq plot of the weighted residuals. A summary plot is also printed: in the diagonal, the value of the weigths vs position of the observations for each root; in the upper diagonal residuals vs residuals of two different roots; in the lower diagonal weights vs weights of two different roots.}

%\details{}

%\value{}

%\references{}

\author{Claudio Agostinelli}

%\note{}

\seealso{\link{wle.lm} a function for estimating linear models with normal distribution error and normal kernel.}

\examples{
library(wle)
# You can find this data set in:
# Hawkins, D.M., Bradu, D., and Kass, G.V. (1984). 
# Location of several outliers in multiple regression data using
# elemental sets. Technometrics, 26, 197-208.
#
data(artificial)

result_wle.lm(y.artificial,x.artificial,boot=40,group=6,num.sol=3)

print(result$coeff)

plot.wle.lm(result)
}

\keyword{Weighted Likelihood}

