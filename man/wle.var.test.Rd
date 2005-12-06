\name{wle.var.test}
\alias{wle.var.test}
\title{Weighted F Test to Compare Two Variances}
\description{
  Performs an Weighted F test to compare the variances of two samples from normal populations. The WF-test is based on weighted likelihood.
}
\usage{
wle.var.test(x, y, ratio = 1, alternative = c("two.sided", "less", "greater"),
         conf.level = 0.95, x.root=1, y.root=1)
}
\arguments{
  \item{x, y}{fitted linear model objects (inheriting from class \code{"wle.lm"}) or fitted normal model objects (inheriting from class \code{"wle.normal"}).}
  \item{ratio}{the hypothesized ratio of the population variances of
    \code{x} and \code{y}.}
  \item{alternative}{the alternative hypothesis; must be one of
    \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.  You
    can specify just the initial letter.}
  \item{conf.level}{confidence level for the returned confidence interval.}
  \item{x.root}{the 'x' root used.}
  \item{y.root}{the 'y' root used.}

}
\details{
  The null hypothesis is that the ratio of the variances in the data to which the normal model (\code{\link{wle.normal}}) or linear models (\code{\link{wle.lm}}) \code{x} and \code{y} were fitted, is equal to \code{ratio}.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the WF test statistic.}
  \item{parameter}{the degrees of the freedom of the WF distribtion of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{conf.int}{a confidence interval for the ratio of the population variances.}
  \item{estimate}{the ratio of the sample variances from \code{x} and \code{y}.}
  \item{null.value}{the ratio of population variances under the null.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{the string \code{"WF test to compare two variances"}.}
  \item{data.name}{a character string giving the names of the data.}
}
%\seealso{}

\references{

Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova (in italian).

Agostinelli, C., (2001) Un approccio robusto alla verifica d'ipotesi basato sulla funzione di verosimiglianza pesata - Robust Testing Hypotheses via Weighted Likelihood function, in press \emph{Statistica}, (in italian).

Agostinelli, C., and Markatou, M., (2001) Test of hypotheses based on the Weighted Likelihood Methodology, \emph{Statistica Sinica},  vol. 11, n. 2, 499-514.
}

\author{Claudio Agostinelli}


\examples{

set.seed(2345)

x <- rnorm(50,0,1)
y <- rnorm(50,10,1)

res.x <- wle.normal(x,group=5)
res.y <- wle.normal(y,group=5)

wle.var.test(res.x, res.y) # Do x and y have the same variance?

set.seed(2345)

x <- c(rnorm(50,0,1),rnorm(20,10,1))
y <- c(rnorm(50,10,1),rnorm(10,0,5))

res.x <- wle.normal(x,group=5,num.sol=2)
res.y <- wle.normal(y,group=5)

res.x
wle.var.test(res.x, res.y, x.root=1)
if (res.x$tot.sol>1) wle.var.test(res.x, res.y, x.root=2)
}

\keyword{robust}
\keyword{htest}


