\name{mle.cv}
\title{Cross Validation Selection Method}
\usage{
mle.cv(formula, data=list(), model=TRUE, x=FALSE, 
       y=FALSE, monte.carlo=500, split, contrasts=NULL)
}
\alias{mle.cv}
\arguments{

  \item{formula}{a symbolic description of the model to be fit.
    The details of model specification are given below.}

  \item{data}{an optional data frame containing the variables
    in the model.  By default the variables are taken from
    the environment which \code{mle.cv} is called from.}

  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding components of the fit (the model frame, the model matrix, the
    response.)}

  \item{monte.carlo}{the number of Monte Carlo replication we use to estimate the average prediction error.}

  \item{split}{the size of the costruction sample. When the suggested value is outside the possible range, the split size is let equal to \eqn{max(round(size^{(3/4)}),nvar+2)}.}

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}
}

\description{The Cross Validation selection method is evaluated for each submodel.}

\details{Models for \code{mle.cv} are specified symbolically.  A typical model has the form \code{response ~ terms} where \code{response} is the (numeric) response vector and \code{terms} is a series of terms which specifies a linear predictor for \code{response}.  A terms specification of the form \code{first+second} indicates all the terms in \code{first} together with all the terms in \code{second} with duplicates removed. A specification of the form \code{first:second} indicates the the set of terms obtained by taking the interactions of all terms in \code{first} with all terms in \code{second}. The specification \code{first*second} indicates the \emph{cross} of \code{first} and \code{second}.  This is the same as \code{first+second+first:second}.
}

\value{ \code{mle.cv} returns an object of \code{\link{class}} \code{"mle.cv"}.

  The function \code{summary} is used to obtain and print a summary of the results.

  The object returned by \code{mle.cv} are:

\item{cv}{the estimated prediction error for each submodels}
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

%\seealso{}

\examples{
library(wle)

data(hald)

cor(hald)

result <- mle.cv(y.hald~x.hald)

summary(result)
}

\keyword{Cross Validation}
\keyworc{CV}
\keyword{model selection}