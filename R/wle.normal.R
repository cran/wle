#############################################################
#                                                           #
#	WLE.NORMAL function                                 #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: December, 19, 2000                            #
#	Version: 0.3                                        #
#                                                           #
#	Copyright (C) 2000 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.normal <- function(x, boot=30, group, num.sol=1, raf="HD", smooth=0.003, tol=10^(-6), equal=10^(-3), max.iter=500)
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

x <- as.vector(x)
size <- length(x)
result <- list()

if (size<2) {
stop("Number of observation must be at least equal to 2")
}

if (group<2) {
group <- max(round(size/4),2)
cat("wle.normal: dimension of the subsample set to default value: ",group,"\n")
}

maxboot <- sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))

if (boot<1 | log(boot) > maxboot) {
stop("Bootstrap replication not in the range")
}

if (!(num.sol>=1)) {
cat("wle.normal: number of solution to report set to 1 \n")
num.sol <- 1
}

if (max.iter<1) {
cat("wle.normal: max number of iteration set to 500 \n")
max.iter <- 500
}

if (smooth<10^(-5)) {
cat("wle.normal: the smooth parameter seems too small \n")
}

if (tol<0) {
cat("wle.normal: the accuracy can not be negative, using default value \n")
tol <- 10^(-6)
}

if (equal<0) {
cat("wle.normal: the equal parameter can not be negative, using default value \n")
equal <- 10^(-3)
}

  z <- .Fortran("wlenorm",
	as.double(x), 
	as.integer(size),
	as.integer(boot),
	as.integer(group),
	as.integer(num.sol),
	as.integer(raf),
	as.double(smooth),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	mean=double(num.sol),
	var=double(num.sol),
	totweight=double(num.sol),
	weight=mat.or.vec(num.sol,size),
	same=integer(num.sol),
	nsol=integer(1),
	nconv=integer(1),
	PACKAGE = "wle")

if (z$nsol>0) {
result$location <- z$mean[1:z$nsol]
result$scale <- sqrt(z$var[1:z$nsol])
     if (z$nsol>1) {
          result$residuals <- matrix(rep(x,z$nsol),nrow=z$nsol,byrow=TRUE) - matrix(rep(z$mean[1:z$nsol],size),nrow=z$nsol,byrow=FALSE)
     } else {
          result$residuals <- x - result$location
     }
result$tot.weights <- z$totweight[1:z$nsol]
result$weights <- z$weight[1:z$nsol,]
result$freq <- z$same[1:z$nsol]
result$tot.sol <- z$nsol
result$not.conv <- z$nconv
result$call <- match.call()

class(result) <- "wle.normal"

return(result)
} else{
stop("No solutions are fuond, checks the parameters")
}
}

print.wle.normal <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Location:\n")
    print.default(format(x$location, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Scale:\n")
    print.default(format(x$scale, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("\nNumber of solutions ",x$tot.sol,"\n")
    cat("\n")
    invisible(x)
}



