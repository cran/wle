#############################################################
#                                                           #
#   wle.wrappednormal function                              #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: Febraury, 3, 2006                                 #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1-4                                           #
#############################################################

wle.wrappednormal <- function(x, mu, rho, sd, K, boot=30, group, num.sol=1, raf="HD", smooth=0.0031, tol=10^(-6), equal=10^(-3), min.sd=1e-3, min.k=10, max.iter=100, use.smooth=TRUE,  p=2, verbose=FALSE) {
  
if (require(circular)) {

    x <- as.circular(x)
    xcircularp <- circularp(x)
    units <- xcircularp$units
    x <- conversion.circular(x, units="radians")
    x <- as.vector(x)
    size <- n <- length(x)
    result <- list()

    sinr <- sum(sin(x))
    cosr <- sum(cos(x))

    est.mu <- FALSE 
    if (missing(mu)) {  
        mu.temp <- atan2(sinr, cosr)
        est.mu <- TRUE
    } else mu.temp <- mu
    
    est.rho <- FALSE
    if (missing(rho)) {
        if (missing(sd)) {
            sd.temp <- sqrt(-2*log(sqrt(sinr^2 + cosr^2)/n))
            if (is.na(sd.temp) || sd.temp < min.sd) sd.temp <- min.sd
            est.rho <- TRUE
         } else sd.temp <- sd
    } else { 
         sd.temp <- sqrt(-2*log(rho))
    }
        
    if (missing(K)) {
        range <- max(mu.temp, x) - min(mu.temp, x)
        K <- (range+6*sd.temp)%/%(2*pi)+1
        K <- max(min.k, K)
    }
    
    
    if (raf!="HD" & raf!="NED" & raf!="SCHI2") stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

    if (missing(group)) {  
        group <- 0
    }

    if (size<2) {
        stop("Number of observation must be at least equal to 2")
    }

    if (group<2) {
        group <- max(round(size/4),1)
        if (verbose) cat("wle.wrappednormal: dimension of the subsample set to default value: ",group,"\n")
    }

    if (boot < 1) {
        stop("Bootstrap replication not in the range")
    }

    if (!(num.sol>=1)) {
        if (verbose) cat("wle.wrappednormal: number of solution to report set to 1 \n")
        num.sol <- 1
    }

    if (max.iter<1) {
        if (verbose) cat("wle.wrappednormal: max number of iteration set to 500 \n")
        max.iter <- 500
    }

    if (tol<=0) {
        if (verbose) cat("wle.wrappednormal: the accuracy must be positive, using default value: 10^(-6) \n")
        tol <- 10^(-6)
    }

    if (equal<=tol) {
        if (verbose) cat("wle.wrappednormal: the equal parameter must be greater than tol, using default value: tol+10^(-3)\n")
        equal <- tol+10^(-3)
    }

    if (p< -1) {
        if (verbose) cat("wle.wrappednormal: the p parameter must be greater than or equal to -1, using default value: 2 \n")
        p <- 2
    }

    if (smooth <= 0) {
        if (verbose) stop("wle.wrappednormal: the smooth parameter must be positive \n")
    }

    tot.sol <- 0
    not.conv <- 0
    iboot <- 0

    while (tot.sol < num.sol & iboot < boot) {
      iboot <- iboot + 1
      continue <- TRUE
      start.iter <- 0
      cat("Bootstrap replication: ", iboot, "\n")
         while (continue & start.iter <= max.iter) {
           x.boot <- sample(x, size=group, replace = TRUE)
#           cat(x.boot, "\n") 
           temp <- mle.wrappednormal(x=x.boot, K=K, tol=tol)
           mu <- temp$mu
           rho <- temp$rho
           sdw <- temp$sd
           continue <- !is.finite(mu) | is.na(mu) | !is.finite(sdw) | is.na(sdw) | sdw <= 0 | !temp$convergence
           start.iter <- start.iter + 1
#           if (verbose) {
#              cat("mu: ", mu, "\n")
#              cat("rho: ", exp(-sdw^2/2), "\n")              
#              cat("sdw: ", sdw, "\n")
#           }
         }
      if (start.iter >= max.iter)
         stop("wle.wrappednormal: the procedure is not able to find a good starting points, please checks the parameters")
     
      if (verbose)
         cat("starting value mu: ", mu, "  rho: ", rho, " sd: ", sdw, "\n")
      x.diff <- tol + 1
      iter <- 0
      while (x.diff > tol & iter < max.iter) {
         iter <- iter + 1 
         mu.old <- mu
         sd.old <- sdw
       
         ff <- density.circular(x=x, z=x, bw=sqrt(smooth)*sdw, adjust=1, kernel="wrappednormal", K=K)$y
  
         if (use.smooth) {
           mm <- dwrappednormal(x, mu=mu, sd=sqrt((1+smooth))*sdw, K=K)
         } else {
           mm <- dwrappednormal(x, mu=mu, sd=sdw, K=K)
         }
         if (any(is.nan(mm)) | any(mm==0)) {
            iter <- max.iter
         } else { 
            dd <- ff/mm - 1

            if (p==Inf && raf=="HD") {
               ww <- log(dd+1)
            } else {
               ww <- switch(raf,
               HD =  p*((dd + 1)^(1/p) - 1) ,
	       NED =  2 - (2 + dd)*exp(-dd) ,
	       SCHI2 =  1-(dd^2/(dd^2 +2)) )       
            }
               
            if (raf=="HD" | raf=="NED") {
               ww <- (ww + 1)/(dd + 1)
            } 

            ww[ww > 1] <- 1
            ww[ww < 0] <- 0
            sommaww <- sum(ww)
                
            iteriter <- 0
            xiterdiff <- 1 + tol
            while (xiterdiff > tol & iteriter <= max.iter) {
              iteriter <- iteriter + 1
              mu.olditer <- mu
              sd.olditer <- sdw

              z <- .Fortran("mlewrpno",
                as.double(x),
                as.double(mu),
                as.double(sdw),
                as.integer(size),
                as.integer(K),
                as.integer(est.mu),
                as.integer(est.rho),
                w=double(size),
                wk=double(size),
                wm=double(size),
                PACKAGE="circular"
              )
              w <- z$w
              wk <- z$wk
              wm <- z$wm
           
              if (est.mu) {
                mu <- c(ww%*%x)/sommaww
                if (any(wk!=0)) {
                  temp <- wk[wk!=0]/w[wk!=0]
                  mu <- mu + 2*pi*c(ww[wk!=0]%*%temp)/sum(ww[wk!=0])
                }
              }
              if (est.rho) {
                if (any(wm!=0)) {
                  temp <- wm[wm!=0]/w[wm!=0]
                  sdw <- sqrt(c(ww[wm!=0]%*%temp)/sum(ww[wm!=0]))
                } else {
                  sdw <- min.sd
                }
              }
              xdiff <- max(abs(mu - mu.olditer), abs(sdw - sd.olditer))
            }
     ########## close the inner while          
           
            x.diff <- max(abs(mu - mu.old), abs(sdw - sd.old))
            if (is.na(x.diff)) iteriter <- iter <- max.iter
          }
     ######## close else for is.na(mm)  
        }
#### end of while (x.diff > tol & iter < max.iter)

#        if (verbose) {
#          cat("mu: ", mu, "\n")
#          cat("rho: ", exp(-sdw^2/2), "\n")              
#          cat("sd: ", sdw, "\n")
#       }
           
   if (iter < max.iter) {
     if (tot.sol==0) {
        mu.store <- mu
        sd.store <- sdw
        w.store <- ww
        tot.store <- sum(ww)/size
        f.store <- ff
        m.store <- mm
        d.store <- dd
        tot.sol <- 1
     } else {
        if (min(abs(mu.store%%(2*pi) - mu%%(2*pi))) > equal | min(abs(sd.store - sdw)) > equal) {
           mu.store <- c(mu.store, mu)
           sd.store <- c(sd.store, sdw)
           w.store <- rbind(w.store, ww)
           tot.store <- c(tot.store, sum(ww)/size)
           f.store <- rbind(f.store, ff)
           m.store <- rbind(m.store, mm)
           d.store <- rbind(d.store, dd)
           tot.sol <- tot.sol + 1
        }
     }
     if (verbose) {
       cat("Number of solutions: ", tot.sol, "\n")
       cat("Last solution found \n")
       cat("mu: ", mu, "\n")
       cat("rho: ", exp(-sdw^2/2), "\n")              
       cat("sd: ", sdw, "\n")
     }
   } else not.conv <- not.conv + 1
 }
##### end of while (tot.sol < num.sol & iboot < boot)
    
 if (tot.sol) {
   result$mu <- c(mu.store)
   result$rho <- c(exp(-sd.store^2/2))
   result$sd <- c(sd.store)
   result$tot.weights <- tot.store
   result$weights <- w.store
   result$f.density <- f.store
   result$m.density <- m.store
   result$delta <- d.store
   result$tot.sol <- tot.sol
   result$not.conv <- not.conv
 } else {
   if (verbose) cat("wle.wrappednormal: No solutions are fuond, checks the parameters\n")
   result$mu <- NA
   result$rho <- NA
   result$sd <- NA
   result$tot.weights <- NA
   result$weights <- rep(NA,size)
   result$f.density <- rep(NA,size)
   result$m.density <- rep(NA,size)
   result$delta <- rep(NA,size)
   result$tot.sol <- 0
   result$not.conv <- boot
 }
     
 result$call <- match.call()
 class(result) <- "wle.wrappednormal"
 return(result)
} else {
  stop("You need package 'circular' for this function")
}
}

#############################################################
#                                                           #
#	print.wle.wrappednormal function                    #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: Febraury, 3, 2006                             #
#	Version: 0.1-1                                      #
#                                                           #
#	Copyright (C) 2006 Claudio Agostinelli              #
#                                                           #
#############################################################

print.wle.wrappednormal <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("mu: ")
    cat(format(x$mu, digits=digits), "\n")
    cat("\n")
    cat("rho: ")    
    cat(format(x$rho, digits=digits), "\n")
    cat("\n")
    cat("sd: ")       
    cat(format(x$sd, digits=digits), "\n")
    cat("\n")   
    cat("\nNumber of solutions ",x$tot.sol,"\n")
    cat("\n")
    invisible(x)
}











