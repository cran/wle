\name{mle.cp}
\title{Mallows Cp}
\usage{
mle.cp(formula, data=list(), model=TRUE, x=FALSE, 
       y=FALSE, var.full=0, contrasts=NULL, verbose=FALSE)
}
\alias{mle.cp}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.
    The details of model specification are given below.}

  \item{data}{an optional data frame containing the variables
    in the model.  By default the variables are taken from
    the environment which \code{mle.cp} is called from.}

  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding components of the fit (the model frame, the model matrix, the
    response).}

  \item{var.full}{the value of variance to be used in the denominator of the Mallows Cp, if 0 the variance estimated from the full model is used.}

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}

  \item{verbose}{if \code{TRUE} warnings are printed.}

}

\description{The Mallows Cp is evaluated for each submodel.}

\details{Models for \code{mle.cp} are specified symbolically.  A typical model has the form \code{response ~ terms} where \code{response} is the (numeric) response vector and \code{terms} is a series of terms which specifies a linear predictor for \code{response}.  A terms specification of the form \code{first+second} indicates all the terms in \code{first} together with all the terms in \code{second} with duplicates removed. A specification of the form \code{first:second} indicates the the set of terms obtained by taking the interactions of all terms in \code{first} with all terms in \code{second}. The specification \code{first*second} indicates the \emph{cross} of \code{first} and \code{second}.  This is the same as \code{first+second+first:second}.
}


\value{ 
    \code{mle.cp} returns an object of \code{\link{class}} \code{"mle.cp"}.

  The function \code{summary} is used to obtain and print a summary of the results, only models below the bisector are reported.
  The generic accessor functions \code{coefficients} and \code{residuals} extract coefficients and residuals returned by \code{mle.cp}.
  The object returned by \code{mle.cp} are:

\item{cp}{Mallows Cp for each submodels}
\item{coefficients}{the parameters estimator, one row vector for eac submodel.}
\item{scale}{an estimation of the error scale, one value for each submodel.}
\item{residuals}{the residuals from the estimated model, one column vector for each submodel.}
\item{call}{the match.call().}
\item{contrasts}{}
\item{xlevels}{}
\item{terms}{the model frame.}
\item{model}{if \code{model=TRUE} a matrix with first column the dependent variable and the remain column the explanatory variables for the full model.}
\item{x}{if \code{x=TRUE} a matrix with the explanatory variables for the full model.}
\item{y}{if \code{y=TRUE} a vector with the dependent variable.}
\item{info}{not well working yet, if 0 no error occurred.}
}

\author{Claudio Agostinelli}

\references{
Mallows, C.L., (1973) Some comments on Cp, \emph{Technometrics}, 15, 661-675.
}

%\seealso{}

\examples{
library(wle)

data(hald)

cor(hald)

result <- mle.cp(y.hald~x.hald)

summary(result)

plot(result)
}

\keyword{regression}


