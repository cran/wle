#############################################################
#                                                           #
#	WLE.NORMAL.MULTI function                           #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: December, 19, 2000                            #
#	Version: 0.3                                        #
#                                                           #
#	Copyright (C) 2000 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.normal.multi <- function(x, boot=30, group, num.sol=1, raf="HD", smooth, tol=10^(-6), equal=10^(-3), max.iter=500)
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

if (is.null(size <- nrow(x)) | is.null(nvar <- ncol(x))) {
    if (is.vector(x)) {
	return(wle.normal(x=x, boot=boot, group=group, num.sol=num.sol, raf=raf, smooth=smooth, tol=tol, equal=equal, max.iter=max.iter)) 
    } else {
	stop("'x' must be a matrix or a vector")
    }
}

if (missing(smooth)) {
    smooth <- wle.smooth(dimension=nvar,costant=4,weight=0.5,interval=c(0.0001,20))$root
}

result <- list()

if (size<(nvar*(nvar+1)/2+nvar)) {
stop(paste("Number of observation must be at least equal to ",nvar*nvar))
}
if (group<(nvar*(nvar+1)/2+nvar)) {
group <- max(round(size/4),(nvar*(nvar+1)/2+nvar))
cat("wle.normal.multi: dimension of the subsample set to default value: ",group,"\n")
}

maxboot <- sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))

if (boot<1 | log(boot) > maxboot) {
stop("Bootstrap replication not in the range")
}

if (!(num.sol>=1)) {
cat("wle.normal.multi: number of solution to report set to 1 \n")
num.sol <- 1
}

if (max.iter<1) {
cat("wle.normal.multi: max number of iteration set to 500 \n")
max.iter <- 500
}

if (smooth<10^(-5)) {
cat("wle.normal.multi: the smooth parameter seems too small \n")
}

if (tol<0) {
cat("wle.normal.multi: the accuracy can not be negative, using default value \n")
tol <- 10^(-6)
}

if (equal<0) {
cat("wle.normal.multi: the equal parameter can not be negative, using default value \n")
equal <- 10^(-3)
}

  z <- .Fortran("wlenormmulti",
	as.double(x), 
	as.integer(size),
	as.integer(nvar),
	as.integer(boot),
	as.integer(group),
	as.integer(num.sol),
	as.integer(raf),
	as.double(smooth),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	mean=mat.or.vec(num.sol,nvar),
	var=mat.or.vec(num.sol,nvar*nvar),
	totweight=double(num.sol),
	weight=mat.or.vec(num.sol,size),
	same=integer(num.sol),
	nsol=integer(1),
	nconv=integer(1))

if (z$nsol>0) {

dn <- colnames(x)

temp <- z$var[1:z$nsol,]
if (z$nsol>1) {
temp.a <- matrix(temp[1,],ncol=nvar)
dimnames(temp.a) <- list(dn,dn)
temp.b <- list(temp.a)
for (i in 2:z$nsol) {
temp.a <- matrix(temp[i,],ncol=nvar)
dimnames(temp.a) <- list(dn,dn)
temp.b <- c(temp.b,list(temp.a))
}
} else {
temp.a <- matrix(temp,ncol=nvar)
dimnames(temp.a) <- list(dn,dn)
temp.b <- list(temp.a)
}

result$location <- z$mean[1:z$nsol,]
result$variance <- temp.b
result$tot.weights <- z$totweight[1:z$nsol]
result$weights <- z$weight[1:z$nsol,]
result$freq <- z$same[1:z$nsol]
result$smooth <- smooth
result$tot.sol <- z$nsol
result$not.conv <- z$nconv
result$call <- match.call()


if (is.null(nrow(result$location))) {
names(result$location) <- dn
} else {
dimnames(result$location) <- list(NULL,dn)
}


class(result) <- "wle.normal.multi"

return(result)
} else {
stop("No solutions are fuond, checks the parameters")
}
}

print.wle.normal.multi <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Location:\n")
    print.default(format(x$location, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("\nVariance-Covariance matrix:\n")
    print.default(x$variance, digits=digits,
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("\nNumber of solutions ",x$tot.sol,"\n")
    cat("\n")
    invisible(x)
}







