\name{wle.t.test}
\alias{wle.t.test}
\alias{print.wle.t.test}
\title{Weighted Likelihood Student's t-Test}
\usage{
wle.t.test(x, y = NULL, alternative = c("two.sided", "less", "greater"),
       mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95,
       boot=30, group, num.sol=1, raf="HD", smooth=0.003, 
       tol=10^(-6), equal=10^(-3), max.iter=500)
}
\arguments{
  \item{x}{a numeric vector of data values.}

  \item{y}{an optional numeric vector data values.}

  \item{alternative}{character specifying the alternative hypothesis,
    must be one of \code{"two.sided"} (default), \code{"greater"} or
    \code{"less"}.  You can specify just the initial letter.}

  \item{mu}{a number indicating the true value of the mean (or
    difference in means if you are performing a two sample test).}

  \item{paired}{a logical indicating whether you want a paired
    weighted t-test.}

  \item{var.equal}{a logical variable indicating whether to treat the
    two variances as being equal. If \code{TRUE} then the pooled
    variance is used to estimate the variance otherwise the Welch
    approximation to the degrees of freedom is used.}

  \item{conf.level}{confidence level of the interval.}

  \item{boot}{the number of starting points based on boostrap subsamples to use in the search of the roots.}

  \item{group}{the dimension of the bootstap subsamples. The default value is \eqn{max(round(size/4),2)} where \eqn{size} is the number of observations.}

  \item{num.sol}{maximum number of roots to be searched.}

  \item{raf}{type of Residual adjustment function to be use:

\code{raf="HD"}: Hellinger Distance RAF,

\code{raf="NED"}: Negative Exponential Disparity RAF,

\code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}

  \item{smooth}{the value of the smoothing parameter.}

  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}

  \item{equal}{the absolute value for which two roots are considered the same. (This parameter must be greater than \code{tol}).}

  \item{max.iter}{maximum number of iterations.}



}
\description{
  wle.t.test performs one and two sample Weighted Likelihood t-tests on vectors of data. This is a robust version of the classical t-test. It should be used when the majority of the data follows a normal distribution.
}
\details{
  If \code{paired} is \code{TRUE} then both \code{x} and \code{y} must
  be specified and they must be the same length.  Missing values are
  removed (in pairs if \code{paired} is \code{TRUE}).  If
  \code{var.equal} is \code{TRUE} then the pooled estimate of the
  variance is used.  By default, if \code{var.equal} is \code{FALSE}
  then the variance is estimated separately for both groups and the
  Welch modification to the degrees of freedom is used.
}
\value{
    The function return a list of class \code{"wle.t.test"} with the following components:

    \item{test}{A list with two dimensions. Each cell is related with a combination of 'x', 'y' roots. In each cell a list of class \code{"htest"} containing the following components:

    \code{statistic} the value of the t-statistic.

    \code{parameters} the degrees of freedom for the t-statistic.

    \code{p.value} the p-value for the test.

    \code{conf.int} a confidence interval for the mean appropriate to the specified alternative hypothesis.

    \code{estimate} the estimated mean or difference in means depending on whether it was a one-sample test or a two-sample test.

    \code{null.value} the specified hypothesized value of the mean or mean difference depending on whether it was a one-sample test or a two-sample test.

    \code{alternative} a character string describing the alternative hypothesis.

    \code{method} a character string indicating what type of t-test was performed.

    \code{data.name} a character string giving the name(s) of the data. 
    
    \code{x.weights} the weights related to the 'x' data.

    \code{y.weights} the weights related to the 'y' data.

    \code{x.root} the number of the 'x' root.

    \code{y.root} the number of the 'y' root.

    }

    \item{x.tot.sol}{the number of solutions for the dataset 'x'.}

    \item{y.tot.sol}{the number of solutions for the dataset 'y' or 1.}

    \item{call}{the match.call().}

    \item{paired}{a logical indicating whether is a paired weighted t-test.}

    \item{x}{'x' data.}

    \item{y}{'y' data or \code{NULL}.}

}

%\seealso{}

\references{

Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D Thesis}, Department of Statistics, University of Padova (in italian).

Agostinelli, C., (1998). Un approccio robusto alla verifica d'ipotesi basato sulla funzione di verosimiglianza pesata - Robust Testing Hypotheses via Weighted Likelihood function, submitted to {\it Statistica}, revised february 2000 (in italian).

Agostinelli, C., Markatou, M. (2001) Test of hypotheses based on the Weighted Likelihood Methodology, in print {\it Statistica Sinica} probably vol. 11, n. 2, Aprile 2001.

}

\author{Claudio Agostinelli}

\examples{
library(wle)
library(ctest)

set.seed(1234)

x <- rnorm(20,0,1)
y <- rnorm(20,6,1)

t.test(x,y)                # P < 2.2e-16
wle.t.test(x,y,group=5)    # P < 2.2e-16

t.test(x,y=c(y,250))       # P = 0.1419 -- NOT significant anymore
wle.t.test(x,y=c(y,250),group=5) # P < 2.2e-16 -- still significant 
set.seed(1234)

# three roots for 'x' and three roots for 'y'
# with nine t-test value
res <- wle.t.test(x=c(rnorm(40,0,1),rnorm(40,10,1)),
           y=c(rnorm(40,0,1),rnorm(40,10,1)),
           group=4,num.sol=3,boot=100)

print(res) # print ALL the t-test
print(res,x.root=1,y.root=1)   # print the test associated to the 
                               # x.root=1,y.root=1 

root.1.1 <- res$test[[1]][[1]] # access to the object associated 
                               # to the x.root=1,y.root=1 

names(root.1.1)

set.seed(1234)

# one root and NOT significant t-test
wle.t.test(x=c(rnorm(40,0,1),rnorm(40,10,1)),
           y=c(rnorm(40,0,1),rnorm(40,10,1)),
           group=4,num.sol=3,boot=100,paired=TRUE)

}

\keyword{robust}
\keyword{htest}