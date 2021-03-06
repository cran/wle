\name{mle.stepwise}
\title{Stepwise, Backward and Forward selection methods}
\usage{
mle.stepwise(formula, data=list(), model=TRUE, x=FALSE, 
             y=FALSE, type="Forward", f.in=4.0, f.out=4.0, 
             contransts=NULL, verbose=FALSE)
}
\alias{mle.stepwise}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.
    The details of model specification are given below.}

  \item{data}{an optional data frame containing the variables
    in the model.  By default the variables are taken from
    the environment which \code{mle.stepwise} is called from.}

  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding components of the fit (the model frame, the model matrix, the
    response.)}

  \item{type}{
\code{type="Stepwise"}: the stepwise methods is used,

\code{type="Forward"}: the forward methods is used,

\code{type="Backward"}: the backward method is used.}

  \item{f.in}{the in value}

  \item{f.out}{the out value}

  \item{contransts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}

  \item{verbose}{if \code{TRUE} warnings are printed.}

}

\description{This function performs Stepwise, Forward and Backward model selection.}

\details{Models for \code{mle.stepwise} are specified symbolically.  A typical model has the form \code{response ~ terms} where \code{response} is the (numeric) response vector and \code{terms} is a series of terms which specifies a linear predictor for \code{response}.  A terms specification of the form \code{first+second} indicates all the terms in \code{first} together with all the terms in \code{second} with duplicates removed. A specification of the form \code{first:second} indicates the the set of terms obtained by taking the interactions of all terms in \code{first} with all terms in \code{second}. The specification \code{first*second} indicates the \emph{cross} of \code{first} and \code{second}.  This is the same as \code{first+second+first:second}.
}

\value{ \code{mle.stepwise} returns an object of \code{\link{class}} \code{"mle.stepwise"}.

  The function \code{summary} is used to obtain and print a summary of the results.

  The object returned by \code{mle.stepwise} are:

\item{step}{the selected models}
\item{type}{the type o model selection procedure was used.}
\item{f.in}{the value of f.in used.}
\item{f.out}{the value of f.out used.}
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

Beale, E.M.L., Kendall, M.G., Mann, D.W., (1967) The discarding of variables in multivariate analysis, \emph{Biometrika}, 54, 357-366.

Efroymson, (1960) Multiple regression analysis, in \emph{Mathematical Methods for Digital Computers}, eds. A. Ralston and H.S. Wilf, 191-203, Wiley, New York.

Garside, M.J., (1965) The best sub-set in multiple regression analysis, \emph{Applied Statistics}, 14, 196-200.

Goldberger, A.S, and Jochems, D.B., (1961) Note on stepwise least squares, \emph{Journal of the American Statistical Association}, 56, 105-110.

Goldberger, A.S., (1961) Stepwise least squares: Residual analysis and specification error, \emph{Journal of the American Statistical Association}, 56, 998-1000.

}

\examples{

library(wle)

data(hald)

cor(hald)

result <- mle.stepwise(y.hald~x.hald)

summary(result)

}

\keyword{regression}




