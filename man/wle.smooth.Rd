\name{wle.smooth}
\title{Bandwidth selection for the normal kernel and normal model.}
\usage{
wle.smooth(weight=0.31, costant=3, level=0.2, dimension=1, raf=1, interval=c(1e-05, 0.5), tol=10^-6, maxiter=1000)
}
\alias{wle.smooth}

\arguments{
 \item{weight}{weights associated to an observation that is \code{costant} scale units away from the mean of the distribution.}
 \item{costant}{Times the contaminated point mass is away from the mean of the distribution in scale units.}
 \item{level}{Mass of the contaminated point.}
 \item{dimension}{Dimension of the normal distribution.}
 \item{raf}{Type of Residual adjustment function to be use:

\code{raf=1}: Hellinger Distance RAF,

\code{raf=2}: Negative Exponential Disparity RAF,

\code{raf=3}: Chi-Squared Disparity RAF.}
 \item{interval}{Interval from which to search the root.}
 \item{tol}{The absolute accuracy to be used to achieve convergence of the algorithm.}
 \item{maxiter}{Maximum number of iterations.}
}
\description{
The bandwidth of the kernel is choose for normal model and normal kernel in such a way a contaminated point \code{costant} times away from the mean of the distribution in scale units and mass \code{level} has a weight no bigger than \code{weight}.
}
\details{
The \code{wle.smooth} use \code{uniroot} function to solve the non linear equation. No handling error is provided yet. For the Chi-Squared Disparity RAFyou should use \code{weight=0.2} and \code{interavl=c(0.1,1)} to have a solution.  
}
\value{A list with four components: root and f.root give the location of the root and the value of the function evaluated at that
point. iter and estim.prec give the number of iterations used and an approximate estimated precision for root.

\code{root} is the value of the bandwidth.

}
\references
{
Agostinelli, C., (1998). Inferenza statistica robusta basata sulla funzione di verosimiglianza pesata: alcuni sviluppi, \emph{ Ph.D. thesis},  Department of Statistics, University of Padova.

Markatou, M., Basu, A. and Lindsay, B.G. (1998). Weighted likelihood estimating equations with a bootstrap root search. \emph{ Journal of the American Statistical Association}, 93, 740-750.

Agostinelli, C., and Markatou, M., (2000). Test of hypotheses based on the Weighted Likelihood Methodology, submitted to \emph{ Statistica Sinica}.
}
\author{Claudio Agostinelli}
%\note{}

\seealso{\code{\link{uniroot}, uniroot}: one dimensional root finding.}

\examples{
wle.smooth()
}

\keyword{Weighted Likelihood}