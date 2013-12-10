#############################################################
#
#	wle.normal.high function
#	Author: Claudio Agostinelli
#	E-mail: claudio@unive.it
#	Date: April, 13, 2011
#	Version: 0.1
#
#	Copyright (C) 2011 Claudio Agostinelli
#
#############################################################

wle.normal.high.old2 <- function(x, boot=30, group, num.sol=1, raf="HD", smooth, tol=10^(-6), equal=10^(-3), max.iter=500, use.smooth=TRUE, tol.int, verbose=FALSE) {

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
    smooth <- 0.008
}

if (!is.logical(use.smooth)) {
    if (verbose) cat("wle.normal.high: the use.smooth must be a logical value, using default value \n")
    use.smooth <- TRUE
}

result <- list()

#if (size<(nvar*(nvar+1)/2+nvar)) {
#stop(paste("Number of observations must be at least equal to ",(nvar*(nvar+1)/2+nvar)))
#}
#if (group<(nvar*(nvar+1)/2+nvar)+1) {
#    group <- max(round(size/4),(nvar*(nvar+1)/2+nvar)+1)
#    if (verbose) cat("wle.normal.high: dimension of the subsample set to default value: ",group,"\n")
#}

maxboot <- sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))

if (boot<1 | log(boot) > maxboot) {
    stop("Bootstrap replication not in the range")
}

if (!(num.sol>=1)) {
    if (verbose) cat("wle.normal.high: number of solution to report set to 1 \n")
    num.sol <- 1
}

if (max.iter<1) {
    if (verbose) cat("wle.normal.high: max number of iteration set to 500 \n")
    max.iter <- 500
}

if (smooth<10^(-5)) {
    if (verbose) cat("wle.normal.high: the smooth parameter seems too small \n")
}

if (tol<=0) {
    if (verbose) cat("wle.normal.high: the accuracy must be positive, using default value: 10^(-6)\n")
    tol <- 10^(-6)
}

if (missing(tol.int)) {
   tol.int <- tol*10^(-4)
} else {
   if (tol.int <=0) {
       if (verbose) cat("wle.normal.high: tol.int must be positive, using default value \n")
       tol.int <- tol*10^(-4) 
   } 
}

if (equal<=tol) {
    if (verbose) cat("wle.normal.high: the equal parameter must be greater than tol, using default value: tol+10^(-3) \n")
    equal <- tol+10^(-3)
}

  result <- list()
  result$freq <- 0
  result$tot.sol <- 0
  result$not.conv <- 0
  iboot <- 0
#####  dsup <- max(x)+12*smooth*nvar
  while (iboot < boot & result$tot.sol < num.sol) {
    iboot <- iboot + 1
    diffp <- diffv <- 1+tol
    xboot <- x[sample(1:size, size=group),]
    posizione <- apply(xboot, 2, median)
    varianza <- diag(apply(xboot, 2, mad))
    iiter <- 0   
    while((diffp > tol | diffv > tol) & iiter < max.iter) {
      iiter <- iiter + 1
      mah <- mahalanobis(x, posizione, varianza)
      posizionevecchia <- posizione
      varianzavecchia <- varianza
      z <- .Fortran("wlechisq",
	as.double(mah),
	as.double(mah),
	as.integer(size),
	as.integer(size),
	as.integer(raf),
        as.double(1),
	as.double(smooth),
        as.integer(1*use.smooth),
        as.double(0),        
        as.double(tol.int),
        as.double(nvar/2),
	weight=double(size),
	density=double(size),
	model=double(size),
	PACKAGE="wle")
      if(iiter == 2) plot(z$weight)
      posizione <- apply(x, 2, function(x) weighted.mean(x, w=z$weight))
      varianza <- cov.wt(x, wt=z$weight)$cov      
      diffp <- max(abs(posizione-posizionevecchia))
      diffv <- max(abs(varianza-varianzavecchia))
    }
    if (iiter < max.iter) {
      if (result$tot.sol==0) {
        result$location <- c(posizione)
        result$variance <- list(varianza)
        result$tot.weights <- sum(z$weight)/size
        result$weights <- z$weight
        result$f.density <- z$density
        result$m.density <- z$model
        result$delta <- z$density/z$model-1
        result$mahalanobis <- mah
        result$freq <- 1
        result$tot.sol <- 1
      } else {
        diffp <- max(abs(t(result$location) - posizione))
        diffv <- rep(0, result$tot.sol)
        for (i in 1:result$tot.sol) {
          diffv[i] <- max(abs(result$variance[[i]]-varianza))
        }
        diffv <- max(diffv)
        if (diffp > equal | diffv > equal) {
          result$location <- rbind(result$location, c(posizione))
          result$variance <- c(result$variance, varianza)
          result$tot.weights <- c(result$tot.weights, sum(z$weight)/size)
          result$weights <- rbind(result$weights, z$weight)
          result$f.density <- rbind(result$f.density, z$density)
          result$m.density <- rbind(result$m.density, z$model)
          result$delta <- rbind(result$delta, z$density/z$model-1)
          result$mahalanobis <- rbind(result$mahalanobis, mah)
          result$freq <- result$freq + 1
          result$tot.sol <- result$tot.sol + 1
        }
      }
    } else {
      result$not.conv <- result$not.conv + 1
    }
  }

dn <- colnames(x)


  if (result$tot.sol==0) {
    if (verbose) cat("wle.normal.high: No solutions are fuond, checks the parameters\n")

    result$location <- rep(NA,nvar)
    result$variance <- matrix(NA,ncol=nvar,nrow=nvar)
    result$tot.weights <- NA
    result$weights <- rep(NA,size)
    result$f.density <- rep(NA,size)
    result$m.density <-rep(NA,size)
    result$delta <- rep(NA,size)
    result$mahalanobis <- rep(NA,size)    
    result$freq <- NA
    result$tot.sol <- 0
    result$not.conv <- boot
  }

if (is.null(nrow(result$location))) {
    names(result$location) <- dn
} else {
    dimnames(result$location) <- list(NULL,dn)
}
result$call <- match.call()
result$smooth <- smooth

class(result) <- "wle.normal.high"

return(result)
}

print.wle.normal.high <- function(x, digits = max(3, getOption("digits") - 3), ...)
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
