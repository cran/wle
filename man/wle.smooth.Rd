\name{wle.smooth}
\title{Bandwidth selection for the normal kernel and normal model.}
\usage{
wle.smooth(weight=0.31,costant=3,level=0.2,
          dimension=1,raf="HD",interval=c(0.00001,0.5),
          tol=10^-6,max.iter=1000)
}
\alias{wle.smooth}
\alias{print.wle.smooth}

\arguments{
  \item{weight}{weights associated to an observation that is \code{costant} scale units away from the mean of the distribution.}

  \item{costant}{times the contaminated point mass is away from the mean of the distribution in scale units.}

  \item{level}{mass of the contaminated point.}

  \item{dimension}{dimension of the normal distribution.}

  \item{raf}{type of Residual adjustment function to be use:

\code{raf="HD"}: Hellinger Distance RAF,

\code{raf="NED"}: Negative Exponential Disparity RAF,

\code{raf="SCHI2"}: Symmetric Chi-Squared Disparity RAF.}

  \item{interval}{interval from which to search the root.}

  \item{tol}{the absolute accuracy to be used to achieve convergence of the algorithm.}

  \item{maxiter}{maximum number of iterations.}
}
\description{
The bandwidth of the kernel is choose for normal model and normal kernel in such a way a contaminated point \code{costant} times away from the mean of the distribution in scale units and mass \code{level} has a weight no bigger than \code{weight}.
}

\details{
The \code{wle.smooth} use \code{uniroot} function to solve the non linear equation. No handling error is provided yet. For the Symmetric Chi-Squared Disparity RAF you should use \code{weight=0.2} and \code{interavl=c(0.1,1)} to have a solution.  
}

\value{
\code{wle.smooth} returns an object of \code{\link{class}} \code{"wle.smooth"}.

  Only print method is implemented for this class.

  The object returned by \code{wle.smooth} is a list with four components: root and f.root give the location of the root and the value of the function evaluated at that
point. iter and estim.prec give the number of iterations used and an approximate estimated precision for root.

\code{root} is the value of the bandwidth.

}

\references{

Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D. thesis},  Department of Statistics, University of Padova.

Markatou, M., Basu, A. and Lindsay, B.G. (1998). Weighted likelihood estimating equations with a bootstrap root search. \emph{ Journal of the American Statistical Association}, 93, 740-750.

Agostinelli, C., and Markatou, M., (2000). Test of hypotheses based on the Weighted Likelihood Methodology, \emph{Statistica Sinica}, in press.
}

\author{Claudio Agostinelli}

\seealso{\code{\link{uniroot}, uniroot}: one dimensional root finding.}

\examples{
library(wle)

wle.smooth()
}

\keyword{robust}


