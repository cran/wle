#############################################################
#                                                           #
#	WLE.CP function                                     #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: December, 19, 2000                            #
#	Version: 0.3                                        #
#                                                           #
#	Copyright (C) 2000 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.cp <- function(formula, data=list(), model=TRUE, x=FALSE, y=FALSE, boot=30, group, var.full=0, num.sol=1, raf="HD", smooth=0.031, tol=10^(-6), equal=10^(-3), max.iter=500, min.weight=0.5, method="full", alpha=2, contrasts=NULL)
{

raf <- switch(raf,
	HD = 1,
	NED = 2,
	SCHI2 = 3,
	-1)

if (raf==-1) stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

type <- switch(method,
	full = 0,
	reduced = 1,
	-1)

if (type==-1) stop("Please, choose the method: full=wieghts based on full model, reduced=weights based on the actual model")

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
    mf$min.weight <- mf$max.iter <- mf$raf <- NULL
    mf$var.full <- mf$alpha <- mf$contrasts <- NULL
    mf$model <- mf$x <- mf$y <- mf$method <- NULL
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

nrep <- 2^nvar-1

if (size<nvar) {
stop("Number of observations must be at least equal to the number of predictors (including intercept)")
}

if (group<nvar) {
group <- max(round(size/4),nvar)
cat("wle.cp: dimension of the subsample set to default value = ",group,"\n")
}

maxboot <- sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))

if (boot<1 | log(boot) > maxboot) {
stop("Bootstrap replication not in the range")
}


if (!(num.sol>=1)) {
cat("wle.cp: number of solution to report set to 1 \n")
num.sol <- 1
}

if (max.iter<1) {
cat("wle.cp: max number of iteration set to 500 \n")
max.iter <- 500
}

if (smooth<10^(-5)) {
cat("wle.cp: the smooth parameter seems too small \n")
}

if (tol<0) {
cat("wle.cp: the accuracy can not be negative, using default value \n")
tol <- 10^(-6)
}
if (equal<0) {
cat("wle.cp: the equal parameter can not be negative, using default value \n")
equal <- 10^(-3)
}

if (var.full<0) {
cat("wle.cp: the variance of the full model can not be negative, using default value \n")
var.full <- 0
}

if (min.weight<0) {
cat("wle.cp: the minimum sum of the weights can not be negative, using default value \n")
min.weight <- 0.5
}

  z <- .Fortran("wlecp",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(0), 
	as.integer(size),
	as.integer(nvar),
	as.integer(boot),
	as.integer(group),
	as.integer(nrep),
	as.integer(raf),
	as.double(smooth),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	as.double(var.full),
	as.integer(num.sol),
	as.double(min.weight),
	as.integer(type),
	as.double(alpha),
	wcp=mat.or.vec(nrep*num.sol,nvar+1),
        param=mat.or.vec(nrep*num.sol,nvar),
	var=double(nrep*num.sol),
	resid=mat.or.vec(nrep*num.sol,size),
	totweight=double(nrep*num.sol),
	weight=mat.or.vec(nrep*num.sol,size),
	same=integer(nrep*num.sol),
	info=integer(1),
	PACKAGE = "wle")

delnull <- z$same==0

result$wcp <- z$wcp[!delnull,]
result$coefficients <- z$param[!delnull,]
result$scale <- sqrt(z$var[!delnull])
result$residuals <- z$resid[!delnull]
result$weights <- z$weight[!delnull,]
result$tot.weights <- z$totweight[!delnull]
result$freq <- z$same[!delnull]
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
dimnames(result$coefficients) <- list(NULL,dn)
dimnames(result$wcp) <- list(NULL,c(dn,"wcp"))

class(result) <- "wle.cp"

return(result)

}

summary.wle.cp <- function (object, num.max=20, ...) {

z <- .Alias(object)
if (is.null(z$terms)) {
    stop("invalid \'wle.cp\' object")
}

if (num.max<1) {
cat("summary.wle.cp: num.max can not less than 1, num.max=1 \n")
num.max <- 1
}

ans <- list()
wcp <- z$wcp
if(is.null(nmodel <- nrow(wcp))) nmodel <- 1
num.max <- min(nmodel,num.max)

if (nmodel!=1) { 
    nvar <- ncol(wcp)-1
    nparam <- apply(wcp[,(1:nvar)],1,sum)
    wcp <- wcp[wcp[,(nvar+1)]<=(nparam+0.00001),]
    if(!is.null(nrow(wcp)) & nrow(wcp)>1) {
	num.max <- min(nrow(wcp),num.max)
    	wcp <- wcp[order(wcp[,(nvar+1)]),]
    	wcp <- wcp[1:num.max,]
    } else num.max <- 1
}

ans$wcp <- wcp
ans$num.max <- num.max
ans$call <- z$call

class(ans) <- "summary.wle.cp"
return(ans)
}

print.wle.cp <- function (x, digits = max(3, getOption("digits") - 3), ...) {
res_summary.wle.cp(object=x, num.max=nrow(x$wcp), ...)
print.summary.wle.cp(res, digits=digits, ...)
}

print.summary.wle.cp <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")

    cat("\nWeighted Mallows Cp:\n")
    if(x$num.max>1) {
    nvar <- ncol(x$wcp)-1
    x$wcp[,(nvar+1)] <- signif(x$wcp[,(nvar+1)],digits)
    } else {
    nvar <- length(x$wcp)-1
    x$wcp[(nvar+1)] <- signif(x$wcp[(nvar+1)],digits)
    }
    print(x$wcp)
    cat("\n")

    cat("Printed the first ",x$num.max," best models \n") 
    invisible(x)
}



