#############################################################
#                                                           #
#	WLE.LM function                                     #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: August, 2, 2001                               #
#	Version: 0.4                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.lm <- function(formula, data=list(), model=TRUE, x=FALSE, y=FALSE, boot=30, group, num.sol=1, raf="HD", smooth=0.031, tol=10^(-6), equal=10^(-3), max.iter=500, contrasts=NULL, verbose=FALSE)
{

raf <- switch(raf,
	HD = 1,
	NED = 2,
	SCHI2 = 3,
	-1)

if (raf==-1) stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

if (missing(group)) {
group <- 0
}

    ret.x <- x
    ret.y <- y
    result <- list()	
    mt <- terms(formula, data = data)
    mf <- cl <- match.call()
    mf$boot <- mf$group <- mf$smooth <- NULL
    mf$tol <- mf$equal <- mf$num.sol <- NULL
    mf$max.iter <- mf$raf <- mf$contrasts <- NULL
    mf$model <- mf$x <- mf$y <- NULL
    mf$verbose <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    xvars <- as.character(attr(mt, "variables"))[-1]
    inter <- attr(mt, "intercept")
    if((yvar <- attr(mt, "response")) > 0) xvars <- xvars[-yvar]
    xlev <-
	if(length(xvars) > 0) {
	    xlev <- lapply(mf[xvars], levels)
	    xlev[!sapply(xlev, is.null)]
	}
    ydata <- model.response(mf, "numeric")
    if (is.empty.model(mt)) 
	stop("The model is empty")
    else 
	xdata <- model.matrix(mt, mf, contrasts)

if (is.null(size <- nrow(xdata)) | is.null(nvar <- ncol(xdata))) stop("'x' must be a matrix")
if (length(ydata)!=size) stop("'y' and 'x' are not compatible")

if (size<nvar) {
    stop("Number of observations must be at least equal to the number of predictors (including intercept)")
}

if (group<nvar) {
    group <- max(round(size/4),nvar)
    if (verbose) cat("wle.lm: dimension of the subsample set to default value: ", group,"\n")
}

maxboot <- sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))

if (boot<1 | log(boot) > maxboot) {
    stop("Bootstrap replication not in the range")
}

if (!(num.sol>=1)) {
    if (verbose) cat("wle.lm: number of solution to report set to 1 \n")
    num.sol <- 1
}

if (max.iter<1) {
    if (verbose) cat("wle.lm: max number of iteration set to 500 \n")
    max.iter <- 500
}

if (smooth<10^(-5)) {
    if (verbose) cat("wle.lm: the smooth parameter seems too small \n")
}

if (tol<=0) {
    if (verbose) cat("wle.lm: the accuracy must be positive, using default value: 10^(-6) \n")
    tol <- 10^(-6)
}

if (equal<=tol) {
    if (verbose) cat("wle.lm: the equal parameter must be greater than tol, using default value: tol+10^(-3) \n")
    equal <- tol+10^(-3)
}

  z <- .Fortran("wleregfix",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(0), 
	as.integer(size),
	as.integer(nvar),
	as.integer(nvar),
	as.integer(boot),
	as.integer(group),
	as.integer(num.sol),
	as.integer(raf),
	as.double(smooth),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	param=mat.or.vec(num.sol,nvar),
	var=double(num.sol),
	resid=mat.or.vec(num.sol,size),
	totweight=double(num.sol),
	weight=mat.or.vec(num.sol,size),
        density=mat.or.vec(num.sol,size),
        model=mat.or.vec(num.sol,size),
        delta=mat.or.vec(num.sol,size),
	same=integer(num.sol),
	nsol=integer(1),
	nconv=integer(1),
	PACKAGE = "wle")

if (z$nsol>0) {
    z$var <- z$var[1:z$nsol]
    z$totweight <- z$totweight[1:z$nsol]
    z$same <- z$same[1:z$nsol]

    if (num.sol==1) {
        z$param <- c(z$param)
        z$resid <- c(z$resid)
        z$weight <- c(z$weight)
        z$density <- c(z$density)
        z$model <- c(z$model)
        z$delta <- c(z$delta)        
    } else {
        z$param <- z$param[1:z$nsol,]
        z$resid <- z$resid[1:z$nsol,]
        z$weight <- z$weight[1:z$nsol,]
        z$density <- z$density[1:z$nsol,]
        z$model <- z$model[1:z$nsol,]
        z$delta <- z$delta[1:z$nsol,]
   }

    y.fit <- t(xdata%*%matrix(z$param,ncol=z$nsol,byrow=TRUE))

    if (z$nsol==1) {
        devparam <- sqrt(z$var*diag(solve(t(xdata)%*%diag(z$weight)%*%xdata,tol=1e-100)))
        y.fit <- as.vector(y.fit)
    } else {
        devparam <- sqrt(z$var[1]*diag(solve(t(xdata)%*%diag(z$weight[1,])%*%xdata,tol=1e-100)))
        for (i in 2:z$nsol) {
             devparam <- rbind(devparam,sqrt(z$var[i]*diag(solve(t(xdata)%*%diag(z$weight[i,])%*%xdata,tol=1e-100))))
        }
    }

    result$coefficients <- z$param
    result$standard.error <- devparam
    result$scale <- sqrt(z$var)
    result$residuals <- z$resid
    result$fitted.values <- y.fit
    result$weights <- z$weight
    result$f.density <- z$density
    result$m.density <- z$model
    result$delta <- z$delta
    result$tot.weights <- z$totweight
    result$tot.sol <- z$nsol
    result$not.conv <- z$nconv
    result$freq <- z$same
} else {
    if (verbose) cat("wle.lm: No solutions are fuond, checks the parameters\n")
    result$coefficients <- rep(NA,nvar)
    result$standard.error <- rep(NA,nvar)
    result$scale <- NA
    result$residuals <- rep(NA,size)
    result$fitted.values <- rep(NA,size)
    result$weights <- rep(NA,size)
    result$f.density <- rep(NA,size)
    result$m.density <- rep(NA,size)
    result$delta <- rep(NA,size)
    result$tot.weights <- NA
    result$tot.sol <- 0
    result$not.conv <- boot
    result$freq <- NA
}

result$call <- cl
result$info <- z$info
result$contrasts <- attr(xdata, "contrasts")
result$xlevels <- xlev
result$terms <- mt

if (model)
    result$model <- mf
if (ret.x)
    result$x <- xdata
if (ret.y)
    result$y <- ydata

dn <- colnames(xdata)

if (is.null(nrow(result$coefficients))) {
names(result$coefficients) <- dn
} else {
dimnames(result$coefficients) <- list(NULL,dn)
}
class(result) <- "wle.lm"

return(result)
}

print.wle.lm <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Coefficients:\n")
    print.default(format(coef(x), digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Scale estimate: ",format(x$scale, digits=digits))
    cat("\n")
    cat("\nNumber of solutions ",x$tot.sol,"\n")
    cat("\n")
    invisible(x)
}

summary.wle.lm <- function(object, root="ALL", ...)
{
z <- .Alias(object)
if (is.null(z$terms)) stop("invalid \'wle.lm\' object")

tot.sol <- z$tot.sol

if (root!="ALL" & !is.numeric(root)) {
    stop("Please, choose one root, for print all root ALL")
} else if (root=="ALL") {
    root <- 1:tot.sol
} else if (tot.sol<root) {
    stop(paste("Root ",root," not found"))
}

ans <- list()
for (iroot in root) {
ans <- c(ans,list(summary.wle.lm.root(object=object, root=iroot)))
}
class(ans) <- "summary.wle.lm" 

return(ans)
}

print.summary.wle.lm <- function (x, digits = max(3, getOption("digits") - 3), signif.stars= getOption("show.signif.stars"),	...)
{
for (i in 1:length(x)) {
print.summary.wle.lm.root(x[[i]], digits = max(3, getOption("digits") - 3), signif.stars= getOption("show.signif.stars"))
}
}




summary.wle.lm.root <- function (object, root=1, ...) {
z <- .Alias(object)
if (is.null(z$terms)) stop("invalid \'wle.lm\' object")

tot.sol <- z$tot.sol
if (tot.sol<root) {
    stop(paste("Root ",root," not found"))
}
if (tot.sol!=1) {
n <- ncol(z$residuals)
p <- ncol(z$coefficients)
} else {
n <- length(z$residuals)
p <- length(z$coefficients)
}

rdf <- z$tot.weights[root]*n - p 
if (tot.sol>1) {   
    r <- z$residuals[root,]
    f <- z$fitted.values[root,]
    w <- z$weights[root,]
    est <- z$coefficients[root,]
    se <- z$standard.error[root,]
} else {
    r <- z$residuals
    f <- z$fitted.values
    w <- z$weights
    est <- z$coefficients
    se <- z$standard.error
}

    mss <- if (attr(z$terms, "intercept")) {
    		m <- sum(w * f /sum(w))
        	sum(w * (f - m)^2)
        	} else sum(w * f^2)
    rss <- sum(w * r^2)
  
    resvar <- rss/rdf
    tval <- est/se
    ans <- z[c("call", "terms")]
    ans$residuals <- sqrt(w)*r
    ans$coefficients <- cbind(est, se, tval, 2*(1 - pt(abs(tval), rdf)))
    dimnames(ans$coefficients)<-
	list(names(est),
	     c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, p)
    if (p != attr(z$terms, "intercept")) {
	df.int <- if (attr(z$terms, "intercept")) 1 else 0
	ans$r.squared <- mss/(mss + rss)
	ans$adj.r.squared <- 1 - (1 - ans$r.squared) *
	    ((n - df.int)/rdf)
	ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
			    numdf = p - df.int, dendf = rdf)
    }

    ans$root <- root
    return(ans)
}

print.summary.wle.lm.root <- function (x, digits = max(3, getOption("digits") - 3), signif.stars= getOption("show.signif.stars"),	...)
{
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
    resid <- x$residuals
    df <- x$df
    rdf <- df[2]
 
   cat("Root ",x$root)
   cat("\n\nWeighted Residuals:\n", sep="")
    if (rdf > 5) {
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	rq <- if (length(dim(resid)) == 2)
	    structure(apply(t(resid), 1, quantile),
		      dimnames = list(nam, dimnames(resid)[[2]]))
	else  structure(quantile(resid), names = nam)
	print(rq, digits = digits, ...)
    }
    else if (rdf > 0) {
	print(resid, digits = digits, ...)
    } else { # rdf == 0 : perfect fit!
	cat("ALL", df[1], "residuals are 0: no residual degrees of freedom!\n")
    }
    if (nsingular <- df[3] - df[1])
	cat("\nCoefficients: (", nsingular,
	    " not defined because of singularities)\n", sep = "")
    else cat("\nCoefficients:\n")

    print.coefmat(x$coef, digits=digits, signif.stars=signif.stars, ...)
    ##
    cat("\nResidual standard error:",
	format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")
    if (!is.null(x$fstatistic)) {
	cat("Multiple R-Squared:", formatC(x$r.squared, digits=digits))
	cat(",\tAdjusted R-squared:",formatC(x$adj.r.squared,d=digits),
	    "\nF-statistic:", formatC(x$fstatistic[1], digits=digits),
	    "on", x$fstatistic[2], "and",
	    x$fstatistic[3], "degrees of freedom,\tp-value:",
	    formatC(1 - pf(x$fstatistic[1], x$fstatistic[2],
			   x$fstatistic[3]), dig=digits),
	    "\n")
    }

    cat("\n")
    invisible(x)
}

fitted.wle.lm <- function(object, ...) object$fitted.values
coef.wle.lm <- function(object, ...) object$coefficients
weights.wle.lm <- .Alias(weights.default)
formula.wle.lm <- function(x, ...) formula(x$terms)

model.frame.wle.lm <-
    function(formula, data, na.action, ...) {
	if (is.null(formula$model)) {
	    fcall <- formula$call
	    fcall$method <- "model.frame"
	    fcall[[1]] <- as.name("wle.lm")
	    eval(fcall, sys.frame(sys.parent()))
	}
	else formula$model
    }

