#############################################################
#                                                           #
#	WLE.GAMMA function                                  #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: February, 16, 2001                            #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.gamma <- function(x, boot=30, group, num.sol=1, raf="HD", smooth=0.008, tol=10^(-6), equal=10^(-3), max.iter=500, shape.int=c(0.01,100), shape.tol=10, use.smooth=TRUE, tol.int) {

wsolve <- function (o, media, medialog) {
   medialog + log(o/media) - digamma(o)
}

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
cat("wle.gamma: dimension of the subsample set to default value: ",group,"\n")
}

maxboot <- sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))

if (boot<1 | log(boot) > maxboot) {
stop("Bootstrap replication not in the range")
}

if (!(num.sol>=1)) {
cat("wle.gamma: number of solution to report set to 1 \n")
num.sol <- 1
}

if (max.iter<1) {
cat("wle.gamma: max number of iteration set to 500 \n")
max.iter <- 500
}

if (smooth<10^(-5)) {
cat("wle.gamma: the smooth parameter seems too small \n")
}

if (tol<0) {
cat("wle.gamma: the accuracy can not be negative, using default value \n")
tol <- 10^(-6)
}

if (equal<0) {
cat("wle.gamma: the equal parameter can not be negative, using default value \n")
equal <- 10^(-3)
}

if (!is.logical(use.smooth)) {
cat("wle.gamma: the use.smooth must be a logical value, using default value \n")
use.smooth <- TRUE
}

if (length(shape.int)!=2) stop("shape.int must be a vector of length 2 \n")

shape.int <- rev(sort(shape.int))
shape.first <- shape.int

if (shape.int[2] <= 0) {
stop("the elements of shape.int must be positive \n")
}

if (shape.int[1] <= 0) {
cat("wle.gamma: the elements of shape.int must be positive, using default value \n")
shape.int[1] <- tol
}

if (shape.tol <=0) {
cat("wle.gamma: shape.tol must be positive, using default value \n")
shape.tol <- 10
} 

if (missing(tol.int)) {
   tol.int <- tol*10^(-4)
} else {
   if (tol.int <=0) {
     cat("wle.gamma: tol.int must be positive, using default value \n")
     tol.int <- tol*10^(-4) 
   } 
}

tot.sol <- 0
not.conv <- 0
iboot <- 0

xlog <- log(x)

while (tot.sol < num.sol & iboot < boot) {
   iboot <- iboot + 1
   x.boot <- x[round(runif(group,0.501,size+0.499))]
   xlog.boot <- log(x.boot)

   media <- sum(x.boot)/group
   medialog<- sum(xlog.boot)/group

   o <- uniroot(wsolve,interval=shape.first, media=media, medialog=medialog)$root

   if (o < tol) o <- 2*tol

   l <- media/o

xdiff <- tol + 1
iter <- 0
while (xdiff > tol & iter < max.iter) {

iter <- iter + 1
shape <- o
lambda <- 1/l
temp <- shape/lambda^2
dsup <- max(x)+ 3*smooth*temp

   z <- .Fortran("wlegamma",
	as.double(x), 
	as.integer(size),
	as.integer(raf),
	as.double(smooth*temp),
        as.integer(1*use.smooth),
        as.double(dsup),
	as.double(tol),
        as.double(tol.int),
	as.double(lambda),
	as.double(shape),
	weights=double(size),
	density=double(size),
	model=double(size))

ww <- z$weights
wsum <- sum(ww)
wmedia <- ww%*%x/wsum
wmedialog <- ww%*%xlog/wsum

shape.int <- c(max(tol,(o-shape.tol)),(o+shape.tol))

o <- uniroot(wsolve, interval=shape.int, media=wmedia, medialog=wmedialog)$root

if (o < tol) o <- 2*tol

l <- wmedia/o

xdiff <- max(abs(c(o-shape,l-1/lambda)))
}

   if (iter < max.iter) {

   if (tot.sol==0) {
      o.store <- o
      l.store <- l
      w.store <- ww
      tot.sol <- 1
   } else {
      if (min(abs(o.store-o))>equal & min(abs(l.store-l))>equal) {
          o.store <- c(o.store,o)
          l.store <- c(l.store,l)
          w.store <- rbind(w.store,ww)
          tot.sol <- tot.sol + 1
      }
   }

   } else not.conv <- not.conv + 1
   

}
##### end of while (tot.sol < num.sol & iboot < boot)



if (tot.sol) {
   result$scale <- c(l.store)
   result$shape <- o.store
   if (tot.sol>1) {
       tot.w <- apply(w.store,1,sum)/size
   } else tot.w <- sum(w.store)/size
  
   result$tot.weights <- tot.w
   result$weights <- w.store
   result$tot.sol <- tot.sol
   result$not.conv <- not.conv
   result$call <- match.call()

   class(result) <- "wle.gamma"

   return(result)
} else{
stop("No solutions are fuond, checks the parameters")
}
}


print.wle.gamma <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Scale:\n")
    print.default(format(x$scale, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Shape:\n")
    print.default(format(x$shape, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("\nNumber of solutions ",x$tot.sol,"\n")
    cat("\n")
    invisible(x)
}






