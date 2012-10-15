\name{residualsAnscombe}
\Rdversion{1.1}
\alias{residualsAnscombe}
\title{
Anscombe residuals
}
\description{
Evaluate the Anscombe residuals for a given type of \code{\link{family}} in GLM.
}

\usage{
residualsAnscombe(y, mu, family, \dots)
}

\arguments{
  \item{y}{vector of the response variable}
  \item{mu}{vector of the same length as \code{y} with the corresponding fitted.values.}
  \item{family}{an object of class \code{\link{family}}.}
  \item{\dots}{not used yet.}
}
\details{
The function performs the Anscombe transformation to obtain residuals that are asymptotically normal distributed. For the Binomial family (see Con and Snell 1968) the transformation is
\deqn{
beta(2/3,2/3)*(pbeta(y/m, 2/3, 2/3) - pbeta(mu-(1-2*mu)/(6*m), 2/3, 2/3))/((mu^(1/6)*(1-mu)^(1/6))/sqrt(m))
}
where \code{m} is the number of trial and \code{y} the number of successes. 
For the Poisson family (see Con and Snell 1968) the transformation is
\deqn{
(3/2*(y^(2/3) - (mu-1/6)^(2/3)))/(mu^(1/6))
}
while for the Gamma family (see McCullagh and Nelder 1989) the transformation is
\deqn{ 
3*(y^(1/3) - mu^(1/3))/(mu^(1/3))
}
and for the Inverse Gaussian family (see McCullagh and Nelder 1989) the transformation is
\deqn{
(\ln(y)-\ln(mu))/\sqrt(mu)
}
}
\value{
It return a vector with the Anscombe residuals.
}
\references{
  Agostinelli, C. and Al-quallaf, F. (2009) Robust inference in 
  Generalized Linear Models. Manuscript in preparation.

  D. R. Cox and E. J. Snell. A general definition of residuals. 
  \emph{Journal of the Royal Statistical Society. Series B} 
  (Methodological), 30(2):248-275, 1968. 

  R. M. Loynes. On cox and snell's general definition of 
  residuals. \emph{Journal of the Royal Statistical Society. Series B} 
  (Methodological), 31(1):103-106, 1969. 

  D. A. Pierce and D. W. Schafer. Residuals in 
  generalized linear models. \emph{Journal of the American Statistical 
  Association}, 81(396):977-986, 1986. 

  Rollin Brant. Residual components in generalized linear 
  models. \emph{The Canadian Journal of Statistics}, 15(2):115-126, 1987. 
}

\author{
Claudio Agostinelli and Fatemah Al-quallaf
}
\seealso{
\code{\link{wle.glm}}
}
\examples{
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
wle.glm.D93 <- wle.glm(counts ~ outcome + treatment, family=poisson())
res <- residualsAnscombe(counts, mu=wle.glm.D93$root1$fitted.values, family=poisson())
qqnorm(res)
qqline(res)
}

\keyword{regression}
