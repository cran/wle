\name{plot.wle.lm}
\title{Plots for the Linear Model}
\usage{
plot.wle.lm(x, level.weight=0.5, 
   ask = interactive() && .Device != "postscript", verbose=FALSE, ...)
}
\alias{plot.wle.lm}

\arguments{
 \item{x}{an object of class \code{wle.lm}}
 \item{level.weight}{value of the weight under which an observations is marked with different color.}
 \item{ask}{logical; if \code{TRUE}, the user is \emph{ask}ed before each plot, see \code{\link{par}(ask=.)}.}
 \item{verbose}{if \code{TRUE} warnings are printed.}
 \item{...}{graphical parameters can be given as arguments.}
}
\description{The \code{plot.wle.lm} function plot a separate graph windows for each roots. In each windows four plots are printed: residuals vs fitted, normal qq plot of the residuals, weighted residuals vs weighted fitted, normal qq plot of the weighted residuals. A summary plot is also printed: in the diagonal, the value of the weigths vs position of the observations for each root; in the upper diagonal residuals vs residuals of two different roots; in the lower diagonal weights vs weights of two different roots.}

%\details{}

%\value{}

%\references{}

\author{Claudio Agostinelli}

%\note{}

\seealso{\code{\link{wle.lm}} a function for estimating linear models with normal distribution error and normal kernel.}

\examples{
library(wle)

data(artificial)

result <- wle.lm(y.artificial~x.artificial,boot=40,group=6,num.sol=3)

result

plot.wle.lm(result)
}

\keyword{regression}
\keyword{robust}

