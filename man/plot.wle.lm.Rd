\name{plot.wle.lm}
\title{Plots for the Linear Model}
\usage{
\method{plot}{wle.lm}(x, roots, which=1:4, which.main, level.weight=0.5, 
   ask = dev.interactive(), col=c(2, 1, 3), id.n=3, labels.id, cex.id = 0.75, verbose=FALSE, ...)
}
\alias{plot.wle.lm}

\arguments{
 \item{x}{an object of class \code{wle.lm}.}
 \item{roots}{a vector specify for which roots the plots are required.}
 \item{which}{if a subset of the plots for each root is required, specify a subset of the numbers \code{0:4}, 0 means no plots.}
 \item{which.main}{if a subset of the plots for the main graphic is
   required, specify a subset of the numbers \code{0:roots^2}, 0 means
   no plots. The plots are specified by columns.}
 \item{level.weight}{value of the weight under which an observations is marked with different color.}
 \item{ask}{logical; if \code{TRUE}, the user is \emph{ask}ed before each plot, see \code{\link{par}(ask=.)}.}
 \item{col}{a vector of 3 elements, to specify colors for the plots.}
 \item{id.n}{number of points to be labelled in some plots, starting
    with the ones with less weight.}
 \item{labels.id}{vector of labels, from which the labels for less weighted
    points will be chosen. If missing uses observation numbers.}
 \item{cex.id}{magnification of point labels.}
 \item{verbose}{if \code{TRUE} warnings are printed.}
 \item{...}{graphical parameters can be given as arguments.}
}
\description{The \code{plot.wle.lm} function plots a separate graph
  windows for each root. In each windows four plots are printed:
  residuals vs fitted, normal qq plot of the residuals, weighted
  residuals vs weighted fitted, normal qq plot of the weighted
  residuals. A summary plot is also printed: in the diagonal, the value
  of the weigths vs position of the observations for each root; in the
  upper diagonal residuals vs residuals of two different roots; in the
  lower diagonal weights vs weights of two different roots. The roots
  and the graphs can be chosen by the arguments \code{roots},
  \code{which.main} and \code{which}.}

\author{Claudio Agostinelli}

\seealso{\code{\link{wle.lm}} a function for estimating linear models with normal distribution error and normal kernel.}

\examples{
library(wle)

data(artificial)

result <- wle.lm(y~x1+x2+x3, data=artificial, boot=40, group=6, num.sol=2)

result

plot.wle.lm(result) # all plots, default behavior 

plot.wle.lm(result, roots=1) # only first root, one plot for window 

par(mfcol=c(2,2))
plot.wle.lm(result, roots=1) # only first root, as usual

plot.wle.lm(result, roots=2, which=1, which.main=0)
# only second root, only residual vs fitted values plot

plot.wle.lm(result, which=1)
# main plot + residual vs fitted values plot for each root

par(mfcol=c(3,2))
plot.wle.lm(result, which=1)
# main plot + residual vs fitted values plot for each root all in the same window


}

\keyword{regression}
\keyword{robust}

