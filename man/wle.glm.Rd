\name{wle.glm}
\Rdversion{1.1}
\alias{wle.glm}
\alias{wle.glm.fit}
\alias{weights.wle.glm}
\alias{print.wle.glm}
\concept{regression}
\concept{logistic}
\concept{log-linear}
\concept{loglinear}


\title{Robust Fitting Generalized Linear Models using Weighted Likelihood}
\description{
\code{wle.glm} is used to robustly fit generalized linear models, specified by
  giving a symbolic description of the linear predictor and a
  description of the error distribution.
}
\usage{
wle.glm(formula, family = binomial, data, weights, subset, 
     na.action, start = NULL, etastart, mustart, offset, 
     control = list(glm = glm.control(...), wle = wle.glm.control()), 
     model = TRUE, method = "wle.glm.fit", x = FALSE, y = TRUE, 
     contrasts = NULL, dist.method = "euclidean", ...)

wle.glm.fit(x, y, weights = NULL, wle.weights = rep(1, NROW(y)), 
     start = NULL, etastart = NULL,  mustart = NULL, offset = rep(0, NROW(y)),
     family = gaussian(), control = list(glm=glm.control(), 
     wle=wle.glm.control()), dist.method='euclidean',
     intercept = TRUE, dispersion=NULL)

\method{weights}{wle.glm}(object, type = c("prior", "working", "wle"), root="all", \dots)
}

\arguments{
  \item{formula}{an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.  The details of model specification are given
    under \sQuote{Details}.}

  \item{family}{a description of the error distribution and link
    function to be used in the model. This can be a character string
    naming a family function, a family function or the result of a call
    to a family function.  (See \code{\link{family}} for details of
    family functions.)}

  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{wle.glm} is called.}

  \item{weights}{an optional vector of \sQuote{prior weights} to be used
    in the fitting process.  Should be \code{NULL} or a numeric vector.}

  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}

  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \sQuote{factory-fresh}
    default is \code{\link{na.omit}}.  Another possible value is
    \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.}

  \item{start}{starting values for the parameters in the linear predictor.}

  \item{etastart}{starting values for the linear predictor.}

  \item{mustart}{starting values for the vector of means.}

  \item{offset}{this can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.
    This should be \code{NULL} or a numeric vector of length equal to
    the number of cases.  One or more \code{\link{offset}} terms can be
    included in the formula instead or as well, and if more than one is
    specified their sum is used.  See \code{\link{model.offset}}.}

  \item{control}{a list with two components of parameters for 
    controlling the fitting process. The first component (\code{glm}) 
    is set using the function \code{\link{glm.control}} while the 
    second component (\code{wle}) is set using the function 
    \code{\link{wle.glm.control}} and it is used to set the parameters 
    regarding the behaviour of the robust method. See the documentation
    of these functions for details.}

  \item{model}{a logical value indicating whether \emph{model frame}
    should be included as a component of the returned value.}

  \item{method}{the method to be used in fitting the model.
    The default method \code{"wle.glm.fit"} uses iteratively reweighted
    least squares (IWLS). The only current alternative is
    \code{"model.frame"} which returns the model frame and does no fitting.}

  \item{x, y}{For \code{wle.glm}:
    logical values indicating whether the response vector and model
    matrix used in the fitting process should be returned as components
    of the returned value.

    For \code{wle.glm.fit}: \code{x} is a design matrix of dimension 
    \code{n * p}, and \code{y} is a vector of observations of length \code{n}.
  }

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}

  \item{dist.method}{
    distance method passed to \code{\link{dist}} to measure the 
    distance between predictor rows.
  }

  \item{intercept}{logical. Should an intercept be included in the
    \emph{null} model?}

  \item{dispersion}{numeric or NULL. If provided used as starting value.}
  
  \item{object}{an object inheriting from class \code{"wle.glm"}.}

  \item{type}{character, partial matching allowed. Type of weights to
    extract from the fitted model object.}

  \item{root}{character ("all") or a number. For which solutions the 
    weights are reported.
  }

  \item{wle.weights}{
    For \code{wle.glm.fit} these are weights used in the iterative algorithm evaluated at each step by \code{wle.glm.weights}.
  }

  \item{\dots}{
    For \code{wle.glm}: arguments to be passed by default to
    \code{\link{glm.control}}: see argument \code{control}.
    
    For \code{weights}:
    further arguments passed to or from other methods.
  }
}

\details{
  A typical predictor has the form \code{response ~ terms} where
  \code{response} is the (numeric) response vector and \code{terms} is a
  series of terms which specifies a linear predictor for
  \code{response}.  For \code{binomial} and \code{quasibinomial}
  families the response can also be specified as a \code{\link{factor}}
  (when the first level denotes failure and all others success) or as a
  two-column matrix with the columns giving the numbers of successes and
  failures.  A terms specification of the form \code{first + second}
  indicates all the terms in \code{first} together with all the terms in
  \code{second} with any duplicates removed.

  A specification of the form \code{first:second} indicates the the set
  of terms obtained by taking the interactions of all terms in
  \code{first} with all terms in \code{second}.  The specification
  \code{first*second} indicates the \emph{cross} of \code{first} and
  \code{second}.  This is the same as \code{first + second +
  first:second}.

  The terms in the formula will be re-ordered so that main effects come
  first, followed by the interactions, all second-order, all third-order
  and so on: to avoid this pass a \code{terms} object as the formula.

  Non-\code{NULL} \code{weights} can be used to indicate that different
  observations have different dispersions (with the values in
  \code{weights} being inversely proportional to the dispersions); or
  equivalently, when the elements of \code{weights} are positive
  integers \eqn{w_i}, that each response \eqn{y_i} is the mean of
  \eqn{w_i} unit-weight observations.  In case of binomial GLM prior weights
  CAN NOT be used to give the number of trials when the response is the
  proportion of successes; in this situation please submit the response 
  variable as two columns (first column successes, second column unsuccesses).
  They would rarely be used for a Poisson GLM.

  \code{wle.glm.fit} is the workhorse function: it is not normally 
  called directly but can be more efficient where the response vector 
  and design matrix have already been calculated. However, this function
  needs starting values and does not look for possible multiple roots in  the system of equations. 

  If more than one of \code{etastart}, \code{start} and \code{mustart}
  is specified, the first in the list will be used.  It is often
  advisable to supply starting values for a \code{\link{quasi}} family,
  and also for families with unusual links such as \code{gaussian("log")}.

  All of \code{weights}, \code{subset}, \code{offset}, \code{etastart}
  and \code{mustart} are evaluated in the same way as variables in
  \code{formula}, that is first in \code{data} and then in the
  environment of \code{formula}.

  For the background to warning messages about \sQuote{fitted probabilities
  numerically 0 or 1 occurred} for binomial GLMs, see Venables &
  Ripley (2002, pp. 197--8).

  Multiple roots may occur if the asymptotic weights are used or in the 
  case of continuous models. The function implements the bootstrap root 
  serach approach described in Markatou, Basu and Lindsay (1998) in order
  to find these roots. 
}

\section{Warnings}{Since in a model selection procedure and/or on an ANOVA table the weights of the WLE procedure must be that of the FULL model (and not that of the actual model) statistics on degrees of freedom, deviance and AIC are valid only if this is the FULL model.}

\value{
  \code{wle.glm} returns an object of class inheriting from 
  \code{"wle.glm"}.

  The function \code{\link{summary}} (i.e., \code{\link{summary.wle.glm}}) can
  be used to obtain or print a summary of the results and the function
  \code{\link{anova}} (i.e., \code{\link{anova.wle.glm.root}})
  to produce an analysis of variance table.

  The generic accessor functions \code{\link{coefficients}},
  \code{effects}, \code{fitted.values} and \code{residuals} can be used to
  extract various useful features of the value returned by \code{wle.glm}.

  \code{weights} extracts a vector of weights, one for each case/root in  the fit (after subsetting and \code{na.action}).

  An object of class \code{"wle.glm"} is a (variable length) list 
  containing at least the following components:

  \code{root1} which is a list with the following components:

    \item{coefficients}{a named vector of coefficients}
    \item{residuals}{the \emph{working} residuals, that is the residuals
    in the final iteration of the IWLS fit.  Since cases with zero
    weights are omitted, their working residuals are \code{NA}.}
    \item{fitted.values}{the fitted mean values, obtained by transforming
    the linear predictors by the inverse of the link function.}
    \item{rank}{the numeric rank of the fitted linear model.}
    \item{family}{the \code{\link{family}} object used.}
    \item{linear.predictors}{the linear fit on link scale.}
    \item{deviance}{up to a constant, minus twice the maximized
    log-likelihood.  Where sensible, the constant is chosen so that a
    saturated model has deviance zero.}
    \item{aic}{Akaike's \emph{An Information Criterion}, minus twice the
    maximized log-likelihood plus twice the number of coefficients (so
    assuming that the dispersion is known).}
    \item{null.deviance}{The deviance for the null model, comparable with
    \code{deviance}. The null model will include the offset, and an
    intercept if there is one in the model.  Note that this will be
    incorrect if the link function depends on the data other than
    through the fitted mean: specify a zero offset to force a correct
    calculation.}
    \item{iter}{the number of iterations of IWLS used.}
    \item{weights}{the \emph{working} weights, that is the weights
    in the final iteration of the IWLS fit.}
    \item{prior.weights}{the weights initially supplied, a vector of
    \code{1}s if none were.}
    \item{df.residual}{the residual degrees of freedom.}
    \item{df.null}{the residual degrees of freedom for the null model.}
    \item{y}{if requested (the default) the \code{y} vector
    used. (It is a vector even for a binomial model.)}
    \item{x}{if requested, the model matrix.}
    \item{model}{if requested (the default), the model frame.}
    \item{converged}{logical. Was the IWLS algorithm judged to have converged?}
    \item{boundary}{logical. Is the fitted value on the boundary of the
    attainable values?}
    \item{wle.weights}{final (robust) weights based on the WLE approach.}
    \item{wle.asymptotic}{logicals. If \code{TRUE} asymptotic weight based on Anscombe residual is used for the corresponding observation.}

    In addition, non-empty fits will have components \code{qr}, 
    \code{R}, \code{qraux}, \code{pivot}
    and \code{effects} relating to the final weighted linear fit.

and the following components:

  \item{family}{the \code{\link{family}} object used.}
  \item{call}{the matched call.}
  \item{formula}{the formula supplied.}
  \item{terms}{the \code{\link{terms}} object used.}
  \item{data}{the \code{data argument}.}
  \item{offset}{the offset vector used.}
  \item{control}{the value of the \code{control} argument used.}
  \item{method}{the name of the fitter function used, currently always
    \code{"wle.glm.fit"}.}
  \item{contrasts}{(where relevant) the contrasts used.}
  \item{xlevels}{(where relevant) a record of the levels of the factors
    used in fitting.}
  \item{tot.sol}{the number of solutions found.}
  \item{not.conv}{the number of starting points that does not converge after the \code{max.iter} (defined using \code{wle.glm.control}) iterations are reached.}
  \item{na.action}{(where relevant) information returned by
    \code{\link{model.frame}} on the special handling of \code{NA}s.}

  Objects of class \code{"wle.glm"} are normally of class 
  \code{"wle.glm"}.

% However, care is needed, as
%  extractor functions for class \code{"wle.glm"} such as
%  \code{\link{residuals}} and \code{weights} do \bold{not} just pick out
%  the component of the fit with the same name.

  If a \code{\link{binomial}} \code{wle.glm} model was specified by 
  giving a two-column response, the weights returned by 
  \code{prior.weights} are
  the total numbers of cases (factored by the supplied case weights) and
  the component \code{y} of the result is the proportion of successes.

  In case of multiple roots (i.e. \code{tot.sol} > 1) then objects of the
  same form as \code{root1} are reported with names \code{root2}, 
  \code{root3} and so on until \code{tot.sol}.
}

\references{
  Agostinelli, C. (1998) Inferenza statistica robusta basata sulla 
  funzione di verosimiglianza pesata: alcuni sviluppi, 
  \emph{Ph.D Thesis}, Department of Statistics, University of Padova.

  Agostinelli, C. and Markatou, M., (1998) A one-step robust estimator 
  for regression based on the weighted likelihood reweighting scheme, 
  \emph{Statistics \& Probability Letters}, Vol. 37, n. 4, 341-350.

  Agostinelli, C. and Markatou, M. (2001) Test of hypotheses based on 
  the Weighted Likelihood Methodology, \emph{Statistica Sinica},  
  vol. 11, n. 2, 499-514.

  Agostinelli, C. and Al-quallaf, F. (2009) Robust inference in 
  Generalized Linear Models. Manuscript in preparation.

  Dobson, A. J. (1990)
  \emph{An Introduction to Generalized Linear Models.}
  London: Chapman and Hall.

  Hastie, T. J. and Pregibon, D. (1992)
  \emph{Generalized linear models.}
  Chapter 6 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

  Markatou, M., Basu, A. and Lindsay, B.G. (1998) Weighted likelihood 
  estimating equations with a bootstrap root search. \emph{Journal of 
  the American Statistical Association}, 93:740-750. 

  McCullagh P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models.}
  London: Chapman and Hall.

  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.}
  New York: Springer.
}

\author{
Claudio Agostinelli and Fatemah Al-quallaf
}

\seealso{
  \code{\link{anova.wle.glm.root}}, \code{\link{summary.wle.glm}}, etc. for
  \code{wle.glm} methods,
  and the generic functions \code{\link{anova}}, \code{\link{summary}},
  \code{\link{effects}}, \code{\link{fitted.values}},
  and \code{\link{residuals}}.

  \code{\link{wle.lm}} for robust non-generalized \emph{linear} models
  for \sQuote{general} linear models.
}

\examples{
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
wle.glm.D93 <- wle.glm(counts ~ outcome + treatment, family=poisson(), x=TRUE, y=TRUE)
wle.glm.D93
anova(extractRoot(wle.glm.D93))
summary(wle.glm.D93)

\dontrun{
## Support for gaussian family not provided yet!
## an example with offsets from Venables & Ripley (2002, p.189)
utils::data(anorexia, package="MASS")

anorex.2 <- wle.glm(Postwt ~ Prewt + Treat + offset(Prewt),
                family = gaussian, data = anorexia)
anorex.2
summary(anorex.2)
}

\dontrun{
# Gamma family is not yet implemented!
# A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100),
    lot1 = c(118,58,42,35,27,25,21,19,18),
    lot2 = c(69,35,26,21,18,16,13,12,12))
wlot1 <- wle.glm(lot1 ~ log(u), data=clotting, family=Gamma,
control=list(glm=glm.control(), wle=wle.glm.control(use.asymptotic=1)))
wlot2 <- wle.glm(lot2 ~ log(u), data=clotting, family=Gamma,
control=list(glm=glm.control(), wle=wle.glm.control(use.asymptotic=1)))
wlot1
wlot2
summary(wlot1)
summary(wlot2)
}
}

\keyword{robust}
\keyword{models}
\keyword{regression}

