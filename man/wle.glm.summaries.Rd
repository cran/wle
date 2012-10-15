\name{wle.glm.summaries}
\alias{family.wle.glm}
\alias{residuals.wle.glm}
\title{Accessing Generalized Linear Model Robust Fits}
\usage{
\method{family}{wle.glm}(object, \dots)

\method{residuals}{wle.glm}(object, type = c("deviance", "pearson", "working",
                           "response", "partial"), root="all", \dots)
}
\arguments{
  \item{object}{an object of class \code{glm}, typically the result of
    a call to \code{\link{glm}}.}
  \item{type}{the type of residuals which should be returned.
    The alternatives are: \code{"deviance"} (default), \code{"pearson"},
    \code{"working"}, \code{"response"}, and \code{"partial"}.}
  \item{root}{vector of integer or characters. If "all" the residuals of all roots are reported, otherwise the position of the root should be supplied.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  These functions are all \code{\link{methods}} for class \code{glm} or
  \code{summary.glm} objects.
}
\details{
  The references define the types of residuals: Davison & Snell is a
  good reference for the usages of each.
  
  The partial residuals are a matrix of working residuals, with each
  column formed by omitting a term from the model.

  How \code{residuals} treats cases with missing values in the original
  fit is determined by the \code{na.action} argument of that fit.
  If \code{na.action = na.omit} omitted cases will not appear in the
  residuals, whereas if \code{na.action = na.exclude} they will appear,
  with residual value \code{NA}.  See also \code{\link{naresid}}.

  For fits done with \code{y = FALSE} the response values are computed
  from other components.
}
\seealso{
  \code{\link{wle.glm}}, \code{\link{anova.wle.glm.root}};
  the corresponding \emph{generic} functions, \code{\link{summary.wle.glm}},
  \code{\link{coef}}, \code{\link{deviance}},
  \code{\link{df.residual}},
  \code{\link{effects}}, \code{\link{fitted}},
  \code{\link{residuals}}.
}

\keyword{models}
\keyword{regression}

