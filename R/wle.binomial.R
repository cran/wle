#############################################################
#                                                           #
#	WLE.BINOMIAL function                               #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: February, 13, 2001                            #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.binomial <- function(x, size, boot=30, group, num.sol=1, raf="HD", tol=10^(-6), equal=10^(-3), max.iter=500)
{

result <- list()

if (raf!="HD" & raf!="NED" & raf!="SCHI2") stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

if (missing(group)) {
group <- 0
}

x <- as.vector(x)
nsize <- length(x)
result <- list()

if (nsize<1) {
stop("Number of observation must be at least equal to 1")
}

if (group<1) {
group <- max(round(nsize/4),1)
cat("wle.binomial: dimension of the subsample set to default value: ",group,"\n")
}

maxboot <- sum(log(1:nsize))-(sum(log(1:group))+sum(log(1:(nsize-group))))

if (boot<1 | log(boot) > maxboot) {
stop("Bootstrap replication not in the range")
}

if (!(num.sol>=1)) {
cat("wle.binomial: number of solution to report set to 1 \n")
num.sol <- 1
}

if (max.iter<1) {
cat("wle.binomial: max number of iteration set to 500 \n")
max.iter <- 500
}

if (tol<0) {
cat("wle.binomial: the accuracy can not be negative, using default value \n")
tol <- 10^(-6)
}

if (equal<0) {
cat("wle.binomial: the equal parameter can not be negative, using default value \n")
equal <- 10^(-3)
}

tot.sol <- 0
not.conv <- 0
iboot <- 0

while (tot.sol < num.sol & iboot < boot) {
   iboot <- iboot + 1
   x.boot <- x[round(runif(group,0.501,nsize+0.499))]
   p <- sum(x.boot)/(size*group)

   ff <- rep(0,nsize)
   x.diff <- tol + 1
   iter <- 0
   while (x.diff > tol & iter < max.iter) {
   iter <- iter + 1
   p.old <- p 
       tff <- table(x)/nsize
       nff <- as.numeric(names(tff))
       for (i in 1:nsize) {
           ff[i] <- tff[nff==x[i]] 
       }
       mm <- dbinom(x,size=size,prob=p)
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

       p <- ww%*%x/(sum(ww)*size)

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
result$p <- p.store
result$tot.weights <- sum(ww)/nsize
result$weights <- w.store
result$tot.sol <- tot.sol
result$not.conv <- not.conv
result$call <- match.call()

class(result) <- "wle.binomial"

return(result)
} else{
stop("No solutions are fuond, checks the parameters")
}
}

print.wle.binomial <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("p:\n")
    print.default(format(x$p, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("\nNumber of solutions ",x$tot.sol,"\n")
    cat("\n")
    invisible(x)
}



