#############################################################
#                                                           #
#	WLE.POISSON function                                #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: February, 16, 2001                            #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.poisson <- function(x, boot=30, group, num.sol=1, raf="HD", tol=10^(-6), equal=10^(-3), max.iter=500)
{

result <- list()

if (raf!="HD" & raf!="NED" & raf!="SCHI2") stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

if (missing(group)) {
group <- 0
}

x <- as.vector(x)
size <- length(x)
result <- list()

if (size<1) {
stop("Number of observation must be at least equal to 1")
}

if (group<1) {
group <- max(round(size/4),1)
cat("wle.poisson: dimension of the subsample set to default value: ",group,"\n")
}

maxboot <- sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))

if (boot<1 | log(boot) > maxboot) {
stop("Bootstrap replication not in the range")
}

if (!(num.sol>=1)) {
cat("wle.poisson: number of solution to report set to 1 \n")
num.sol <- 1
}

if (max.iter<1) {
cat("wle.poisson: max number of iteration set to 500 \n")
max.iter <- 500
}

if (tol<0) {
cat("wle.poisson: the accuracy can not be negative, using default value \n")
tol <- 10^(-6)
}

if (equal<0) {
cat("wle.poisson: the equal parameter can not be negative, using default value \n")
equal <- 10^(-3)
}

tot.sol <- 0
not.conv <- 0
iboot <- 0

while (tot.sol < num.sol & iboot < boot) {
   iboot <- iboot + 1
   x.boot <- x[round(runif(group,0.501,size+0.499))]
   p <- sum(x.boot)/group

   ff <- rep(0,size)
   x.diff <- tol + 1
   iter <- 0
   while (x.diff > tol & iter < max.iter) {
   iter <- iter + 1
   p.old <- p 
       tff <- table(x)/size
       nff <- as.numeric(names(tff))
       for (i in 1:size) {
           ff[i] <- tff[nff==x[i]] 
       }
       mm <- dpois(x,lambda=p)
       dd <- ff/mm - 1
       
       ww <- switch(raf,
                 HD =  2*(sqrt(dd + 1) - 1) ,
	         NED =  2 - (2 + dd)*exp(-dd) ,
	         SCHI2 =  1-(dd^2/(dd^2 +2)) )       

       if (raf=="HD" | raf=="NED") {
            ww <- (ww + 1)/(dd + 1)
       }

       ww[ww > 1] <- 1
       ww[ww < 0] <- 0

       p <- ww%*%x/sum(ww)

       x.diff <- abs(p - p.old)
   }
#### end of while (x.diff > tol & iter < max.iter)

   if (iter < max.iter) {

   if (tot.sol==0) {
      p.store <- p
      w.store <- ww
      tot.sol <- 1
   } else {
      if (min(abs(p.store-p))>equal) {
          p.store <- c(p.store,p)
          w.store <- rbind(w.store,ww)
          tot.sol <- tot.sol + 1
      }
   }

   } else not.conv <- not.conv + 1
   

}
##### end of while (tot.sol < num.sol & iboot < boot)

if (tot.sol) {
result$lambda <- c(p.store)
result$tot.weights <- sum(ww)/size
result$weights <- w.store
result$tot.sol <- tot.sol
result$not.conv <- not.conv
result$call <- match.call()

class(result) <- "wle.poisson"

return(result)
} else{
stop("No solutions are fuond, checks the parameters")
}
}

print.wle.poisson <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("lambda:\n")
    print.default(format(x$lambda, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("\nNumber of solutions ",x$tot.sol,"\n")
    cat("\n")
    invisible(x)
}



