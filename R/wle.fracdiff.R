#############################################################
#                                                           #
#	wle.fracdiff.ao function                                #
#	Author: Claudio Agostinelli                             #
#	E-mail: claudio@unive.it                                #
#	Date: April, 08, 2002                                   #
#	Version: 0.2-2                                          #
#                                                           #
#     Copyright (C) 2002 Claudio Agostinelli                #
#                                                           #
#############################################################

wle.fracdiff.ao <- function(d, sigma2, x, M=100, x.init=rep(0,M), x.mean=0, use.init=FALSE, raf=1, smooth=0.0031, w.level=0.5, verbose=FALSE, ao.list=list(0), popolation.size=20, popolation.choose=5, elements.random=4, num.max=length(x)) {

    if (use.init) {
        MM <- 0
    } else {
        MM <- M
    }

    nused <- length(x)
    resid <- wle.fracdiff.residuals(d=d, M=M, x=x, x.ao=x, x.init=x.init, x.mean=x.mean, use.init=use.init)  
    nresid <- length(resid)

    weights <- .Fortran("wlew",
	as.double(resid), 
	as.integer(nresid),
	as.double(resid), 
	as.integer(nresid), 
	as.integer(raf),
	as.double(smooth),
	as.double(sigma2),
	totweight=double(1),
	weights=double(nresid),
	PACKAGE="wle")$weights

    ao.position <- 0
    pos.temp <- 1:nresid
    pos.temp <- pos.temp[order(weights)]
    weights.sort <- sort(weights)
    ao.temp <- weights.sort <= w.level
    pos.temp <- pos.temp[ao.temp]+MM
    ao <- rep(FALSE,nused)
    if (length(pos.temp)) {
        pos.temp <- pos.temp[1:min(length(pos.temp),num.max)]
        ao[pos.temp] <- TRUE
    }

    pos <- (1:nused)[ao]

    if (verbose) {
        cat("We have the following observations under the w.level=",w.level,":\n",pos,"\n")
    }

    if (any(ao)) {
        model.in <- vector(length=0)

        for (i in 1:length(ao.list)) {
             if (all(is.element(ao.list[[i]],pos))) {
                 temp <- vector(length=0)
                 for (j in 1:length(ao.list[[i]])) {
                      temp <- c(temp,(1:length(pos))[pos==ao.list[[i]][j]])
                 }
                 model.in <- c(model.in,sum(2^(temp-1)))
             }
         }

         num.model <- max(length(model.in),popolation.size)
         num.pos <- (2^sum(ao))-1
         dim.dim <- floor(log(num.pos,2))+1
         w.tilde <- rep(0,num.model)

         model.in <- c(model.in,wle.riunif((num.model-length(model.in)),1,num.pos))

         for (isearch in 1:num.model) {
              pos.ao <- sort(pos[binary(model.in[isearch],dim.dim)$dicotomy])
              num.ao <- length(pos.ao)
              x.ao <- x
              for (t in pos.ao) {
                   x.ao[t] <- wle.fracdiff.fitted(t=t, d=d, M=M, x=x.ao, x.init=x.init, x.mean=x.mean, use.init=use.init)
              }
              resid.ao <- wle.fracdiff.residuals(d=d, M=M, x=x, x.ao=x.ao, x.init=x.init, x.mean=x.mean, use.init=use.init) 
              resid.ao <- resid.ao[-pos.ao]
              w.temp <- wle.weights(x=resid.ao, smooth=smooth, sigma2=sigma2, raf=raf, location=TRUE)

              weights.temp <- w.temp$weights
              w.tilde[isearch] <- sum(weights.temp)/nresid
         }

         model.in <- model.in[order(w.tilde)]
         w.tilde <- sort(w.tilde)

         while ((model.in[1]-model.in[num.model])!=0) {
                num.model.sel <- popolation.choose
                cum.wtilde <- cumsum(w.tilde)[num.model.sel:num.model]
                pos.child <- vector(length=0)

                while (length(pos.child)==0) {
                       temp <- runif(2,0,cum.wtilde[length(cum.wtilde)])
                       pos.aaa <- min((num.model.sel:num.model)[cum.wtilde > temp[1]])
                       pos.bbb <- min((num.model.sel:num.model)[cum.wtilde > temp[2]])

                       pos.aa <- pos[binary(model.in[pos.aaa],dim.dim)$dicotomy]
                       pos.bb <- pos[binary(model.in[pos.bbb],dim.dim)$dicotomy]

                       pos.child <- c(pos.aa,pos.bb,pos[wle.riunif(elements.random,1,length(pos))])
                       pos.child <- pos.child[as.logical(wle.riunif(length(pos.child),0,1))]
                       pos.child <- sort(unique(pos.child))
                }

                temp.child <- vector(length=0)
                for (i in 1:length(pos.child)) {
                     temp.child <- c(temp.child,(1:length(pos))[pos==pos.child[i]])
                }

                model.child <- sum(2^(temp.child-1))

                num.child <- length(pos.child)
                x.ao <- x
                for (t in pos.child) {
                     x.ao[t] <- wle.fracdiff.fitted(t=t, d=d, M=M, x=x.ao, x.init=x.init, x.mean=x.mean, use.init=use.init)
                }

                resid.ao <- wle.fracdiff.residuals(d=d, M=M, x=x, x.ao=x.ao, x.init=x.init, x.mean=x.mean, use.init=use.init)
                resid.ao <- resid.ao[-pos.child]

                w.temp <- wle.weights(x=resid.ao, smooth=smooth, sigma2=sigma2, raf=raf, location=TRUE)

                weights.temp <- w.temp$weights
                w.tilde.child <- sum(weights.temp)/nresid

                w.tilde <- c(w.tilde,w.tilde.child)
                model.in <- c(model.in,model.child)

                model.in <- model.in[order(w.tilde)][-1]
                w.tilde <- sort(w.tilde)[-1]
         }

         if (max(w.tilde)<(sum(weights)/nresid)) {
             ao.position <- NULL
         } else {
             ao.position <- sort(pos[binary(model.in[1],dim.dim)$dicotomy])
         }

    } else {
        ao.position <- NULL
    }

    x.ao <- x
    for (t in ao.position) {
         x.ao[t] <- wle.fracdiff.fitted(t=t, d=d, M=M, x=x.ao, x.init=x.init, x.mean=x.mean, use.init=use.init)
    }

    resid.ao <- wle.fracdiff.residuals(d=d, M=M, x=x, x.ao=x.ao, x.init=x.init, x.mean=x.mean, use.init=use.init)
    w.temp <- wle.weights(x=resid.ao, smooth=smooth, sigma2=sigma2, raf=raf, location=TRUE)
    resid.ao <- resid.ao - w.temp$location

    if (verbose) {
        cat("Additive outliers: \n", ao.position, "\n")
    }

    return(x.ao=x.ao, resid.ao=resid.ao, ao.position=ao.position)
}

#############################################################
#                                                           #
#	wle.fracdiff function                              	    #
#	Author: Claudio Agostinelli                             #
#	E-mail: claudio@unive.it                                #
#	Date: April, 02, 2002                                   #
#	Version: 0.1-2                                          #
#                                                           #
#	Copyright (C) 2002 Claudio Agostinelli                  #
#                                                           #
#############################################################

wle.fracdiff <- function(x, lower, upper, M, group, na.action=na.fail, tol=10^(-6), equal=10^(-3), raf="HD", smooth=0.0031, smooth.ao=smooth, boot=10, num.sol=1, x.init=rep(0,M), use.uniroot=FALSE, use.init=FALSE, max.iter.out=20, max.iter.in=100, max.iter.step=5000, max.iter.start=max.iter.step,  verbose=FALSE, w.level=0.4, min.weights=0.5, popolation.size=10, popolation.choose=5, elements.random=2, init.values=NULL, num.max=length(x), include.mean=FALSE, ao.list=list(0)) {

    if (use.init) {
        MM <- 0
    } else {
        MM <- M
    }

    raf <- switch(raf,
	HD = 1,
	NED = 2,
	SCHI2 = 3,
	-1)

    if (raf==-1) stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

    result <- list()
    series <- deparse(substitute(x))
    if(NCOL(x) > 1) stop("only implemented for univariate time series")

    x <- na.action(as.ts(x))
    nused <- length(x)

    if (length(x.init)!=M) stop("x.init must have M elements\n")

# start bootstrap iteration
    first.time <- TRUE
    iboot <- 1
    tot.sol <- 0
    not.conv <- 0

    while (iboot<=boot & tot.sol<num.sol) {
           pos.iboot <- round(runif(1,(group+1),nused))
           x.boot <- x[(pos.iboot-group+1):pos.iboot]

           if (!is.null(init.values)) {
               temp <- list(conv=TRUE)
               temp$d <- init.values[1]
               temp$sigma2 <- init.values[2]
               temp$x.mean <- init.values[3]
               temp$resid <- wle.fracdiff.residuals(d=temp$d, M=M, x=x, x.ao=x, x.init=x.init, x.mean=temp$x.mean, use.init=use.init)
           } else {
               temp <- wle.fracdiff.solve(x=x.boot, x.init=x.init, max.iter=max.iter.start, verbose=verbose, M=M,  lower, upper, tol=tol, use.uniroot=use.uniroot, use.init=use.init, include.mean=include.mean)
               temp$resid <- wle.fracdiff.residuals(temp$d, M=M, x=x, x.ao=x, x.init=x.init, x.mean=temp$x.mean, use.init=use.init)
               temp$sigma2 <- wle.fracdiff.sigma2(resid=temp$resid)
           }
    
           if (temp$conv) {
               d <- temp$d
               sigma2 <- temp$sigma2
               x.mean <- temp$x.mean
               resid <- temp$resid
               nresid <- length(resid)

               if (verbose) {
	           cat("Initial values from the subsample ",iboot,": \n parameters, d: ", d,"\n sigma2: ",sigma2," \n x.mean: ",x.mean, " \n") 
               }

               weights <- .Fortran("wlew",
    	            as.double(resid),
    	            as.integer(nresid),
	            as.double(resid),
	            as.integer(nresid),
	            as.integer(raf),
	            as.double(smooth.ao),
	            as.double(sigma2),
	            totweights=double(1),
	            weights=double(nresid),
				PACKAGE="wle")$weights

               if (sum(weights)/nresid >= min.weights) {
                   wres <- wle.fracdiff.ao(d=d, sigma2=sigma2, x=x, M=M, x.init=x.init, x.mean=x.mean, use.init=use.init, raf=raf, smooth=smooth.ao, w.level=w.level, verbose=verbose, ao.list=ao.list, popolation.size=popolation.size, popolation.choose=popolation.choose, elements.random=elements.random, num.max=num.max)

                   x.ao <- wres$x.ao
                   ao.position <- wres$ao.position
                   resid <- wle.fracdiff.residuals(d, M=M, x=x, x.ao=x.ao, x.init=x.init, x.mean=x.mean, use.init=use.init)
                   weights <- wle.weights(x=resid, smooth=smooth, sigma2=sigma2, raf=raf, location=TRUE)$weights

                   if (!is.null(ao.position)) {
                       ao.list <- c(ao.list,list(ao.position))
                   }
                   ao.position.old <- c(ao.position,0)
                   conv <- TRUE
                   iter.out <- 0

                   while (!setequal(ao.position,ao.position.old) & conv) {
    	                  iter.out <- iter.out + 1
    	                  ao.position.old <- ao.position	
                          maxtol <- tol + 1
                          iter.in <- 0
                          while (maxtol > tol & conv) {
                                 iter.in <- iter.in + 1
                                 d.old <- d
                                 sigma2.old <- sigma2
                                 x.mean.old <- x.mean
	                         res <- wle.fracdiff.solve(x=x.ao, M=M, x.init=x.init, lower=lower, upper=upper, w=weights, tol=tol, max.iter=max.iter.step, verbose=verbose, use.uniroot=use.uniroot, use.init=use.init, include.mean=include.mean)
	                         d <- res$d
                                 x.mean <- res$x.mean
                                 resid <- wle.fracdiff.residuals(d=d, M=M, x=x, x.ao=x.ao, x.init=x.init, x.mean=x.mean, use.init=use.init)
	                         sigma2 <- wle.fracdiff.sigma2(resid=resid, w=weights)
                                 weights <- wle.weights(x=resid, smooth=smooth, sigma2=sigma2, raf=raf, location=TRUE)$weights
                                 conv <- res$conv            

                                 if (iter.in > max.iter.in) {
                                     if (verbose) cat("Convergence problem: maximum iteration number reached in the outer loop\n")
                                     conv <- FALSE
                                 }
                                 maxtol <- max(abs(d-d.old),abs(sigma2-sigma2.old), abs(x.mean-x.mean.old))
                                  if (verbose)  {
	    	                      cat("inner loop, iteration: ",iter.in," \n parameters, d: ",d," \n sigma2: ",sigma2," \n x.mean: ",x.mean," \n")
	                          }
                          }

                          if (conv) {      
	                      if (verbose) {
	    	                  cat("outer loop, iteration: ",iter.out," convergence achieved for the inner loop \n")
	                      }

                              wres <- wle.fracdiff.ao(d=d, sigma2=sigma2, x=x, M=M, x.init=x.init, x.mean=x.mean, use.init=use.init, raf=raf, smooth=smooth.ao, w.level=w.level, verbose=verbose, ao.list=ao.list, popolation.size=popolation.size, popolation.choose=popolation.choose, elements.random=elements.random, num.max=num.max)

                              x.ao <- wres$x.ao
                              ao.position <- wres$ao.position
                              resid <- wle.fracdiff.residuals(d, M=M, x=x, x.ao=x.ao, x.init=x.init, x.mean=x.mean, use.init=use.init)
                              sigma2 <- wle.fracdiff.sigma2(resid=resid, w=weights)
                              weights <- wle.weights(x=resid, smooth=smooth, sigma2=sigma2, raf=raf, location=TRUE)$weights

                              if (!is.null(ao.position)) {
                                  ao.list <- c(ao.list,list(ao.position))
                              }

                              if ((lao <- length(ao.list))>2) {
                                  for (i.ao in 1:(lao-2)) {
                                       if (setequal(ao.list[[lao]],ao.list[[i.ao]])) {
                                           conv <- FALSE
                                       }
                                   } 
                               }

                               if (iter.out > max.iter.out) {
                                   if (verbose) cat("Convergence problem: maximum iteration number reached in the outer loop\n")
                                   conv <- FALSE
                               }

                          }

                   }
# end while (!setequal(ao.position,ao.position.old) & conv)

                   if (conv) {
                       resid.with.ao <- wle.fracdiff.residuals(d, M=M, x=x, x.ao=x, x.init=x.init, x.mean=x.mean, use.init=use.init) 
                       resid.with.ao <- ts(resid.with.ao, start=(start(x)+MM), end=end(x), frequency=frequency(x))
                       class(resid.with.ao) <- "ts"

                       weights.with.ao <- wle.weights(x=resid.with.ao, smooth=smooth, sigma2=sigma2, raf=raf, tol=tol, location=TRUE)$weights

                       resid <- ts(resid, start=(start(x)+MM), end=end(x), frequency=frequency(x))        
                       class(resid) <- "ts" 

                       resid.without.ao <- wle.fracdiff.residuals(d, M=M, x=x.ao, x.ao=x.ao, x.init=x.init, x.mean=x.mean, use.init=use.init) 
                       resid.without.ao <- ts(resid.without.ao, start=(start(x)+MM), end=end(x), frequency=frequency(x))
                       class(resid.without.ao) <- "ts"
 
                       weights.without.ao <- wle.weights(x=resid.without.ao, smooth=smooth, sigma2=sigma2, raf=raf, tol=tol, location=TRUE)$weights
 
                       x.ao <- ts(x.ao, start=start(x), end=end(x), frequency=frequency(x))
                       class(x.ao) <- "ts" 

    if (first.time) {
        d.final <- d
        sigma2.final <- sigma2
        x.mean.final <- c(x.mean)
        weights.final <- weights
        weights.with.ao.final <- weights.with.ao
        weights.without.ao.final <- weights.without.ao   
        resid.final <- resid
        resid.with.ao.final <- resid.with.ao
        resid.without.ao.final <- resid.without.ao
        x.ao.final <- x.ao
        ao.position.final <- list(ao.position)
        first.time <- FALSE	
        tot.sol <- 1
    } else {
        if (min(abs(d.final-d))>equal) {
	    tot.sol <- tot.sol+1
	    d.final <- c(d.final,d)
            sigma2.final <- c(sigma2.final,sigma2)
            x.mean.final <- c(x.mean.final,c(x.mean))
	    weights.final <- rbind(weights.final,weights)
            weights.with.ao.final <- rbind(weights.with.ao.final, weights.with.ao)
            weights.without.ao.final <- rbind(weights.without.ao.final, weights.without.ao)
	    resid.final <- rbind(resid.final,resid)
            resid.with.ao.final <- rbind(resid.with.ao.final,resid.with.ao)
            resid.without.ao.final <- rbind(resid.without.ao.final,resid.without.ao)
            ao.position.final <- c(ao.position.final,list(ao.position))
            x.ao.final <- rbind(x.ao.final,x.ao)
        }
    }

                   } else {
    	               not.conv <- not.conv+1
                   }
               } else {
    	           not.conv <- not.conv+1
               }
           } else {
    	       not.conv <- not.conv+1
           }    

        iboot <- iboot+1
    }
# end bootstrap iteration

if(tot.sol==0) {
   
    result$d <- NULL
    result$sigma2 <- NULL
    result$x.mean <- NULL
    result$resid <- NULL
    result$resid.with.ao <- NULL
    result$resid.without.ao <- NULL
    result$x.ao <- NULL
    result$call <- match.call()
    result$weights <- NULL
    result$weights.with.ao <- NULL
    result$weights.without.ao <- NULL   
    result$tot.sol <- 0
    result$not.conv <- not.conv
    result$ao.position <- NULL
} else { 
    result$d <- d.final
    result$sigma2 <- sigma2.final
    result$x.mean <- x.mean.final
    result$resid <- resid.final
    result$resid.without.ao <- resid.without.ao.final
    result$resid.with.ao <- resid.with.ao.final
    result$x.ao <- x.ao.final
    result$call <- match.call()
    result$weights <- weights.final
    result$weights.with.ao <- weights.with.ao.final
    result$weights.without.ao <- weights.without.ao.final
    result$tot.sol <- tot.sol
    result$not.conv <- not.conv
    result$ao.position <- ao.position.final
    }

    class(result) <- "wle.farima"
    
return(result)
}


#############################################################
#                                                           #
#	wle.fracdiff.solve function                         #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: December, 11, 2001                            #
#	Version: 0.1-1                                      #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fracdiff.solve <- function(x, M=100, x.init=rep(0,M), lower, upper, w=rep(1,length(x)), tol=.Machine$double.eps^0.25, max.iter=1000, verbose=FALSE, use.uniroot=FALSE, use.init=FALSE, include.mean=FALSE) {

    result <- list()

    if (include.mean) {
        x.mean <- w%*%x/sum(w)
        x <- x - x.mean
    } else {
        x.mean <- 0
    }

    result$x.mean <- x.mean
    result$call <- match.call()

    if (use.uniroot) {
        temp <- uniroot(wle.fracdiff.equation, x=x, M=M, w=w, x.init=x.init, lower=lower, upper=upper, use.uniroot=use.uniroot, tol=tol, maxiter=max.iter, verbose=verbose, use.init=use.init)
        if (temp$iter < max.iter) {
            result$d <- temp$root
            result$conv <- TRUE
        } else {
            result$d <- NA
            result$conv <- FALSE
        }
    } else {
        temp <- optimize(wle.fracdiff.equation, x=x, M=M, w=w, x.init=x.init, lower=lower, upper=upper, use.uniroot=use.uniroot, tol=tol, verbose=verbose, use.init=use.init)
        result$d <- temp$minimum
        result$conv <- TRUE
    }
    return(result)
}

#############################################################
#                                                           #
#	wle.fracdiff.equation function                      #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: November, 30, 2001                            #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################


wle.fracdiff.equation <- function(d, M, x, x.init=rep(0,M), w=rep(1,length(x)), use.uniroot=FALSE, use.init=FALSE, verbose=FALSE) {
 
    if (use.init) {
        MM <- 0
    } else {
        MM <- M
    }
    nused <- length(x)
    pi.coef <- wle.fracdiff.pi.coef(d,M)
    if (use.uniroot) {
        xi.coef <- wle.fracdiff.xi.coef(d,M)
    }
    y <- c(x.init,x)

    somma <- 0
    if (use.uniroot) {
        for (t in (1+MM):nused) {
             somma <- somma + w[t]*(x[t]+pi.coef%*%y[(t-1+M):t])*(xi.coef%*%y[(t-1+M):t])
        }
    } else {
        for (t in (1+MM):nused) {
             somma <- somma + w[t]*(x[t]+pi.coef%*%y[(t-1+M):t])^2
        }
    }

    if (verbose) cat("value of d: ",d," value of the function: ",somma,"\n")

    return(as.vector(somma))
}


#############################################################
#                                                           #
#	wle.fracdiff.pi.coef function                       #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: November, 30, 2001                            #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fracdiff.pi.coef <- function(d,M) {
     pi.coef <- rep(0,M)
     pi.coef[1] <- -d
     for (j in 2:M) {
          pi.coef[j] <- pi.coef[j-1]*(j-1-d)/j
     }
return(pi.coef)
}

#############################################################
#                                                           #
#	wle.fracdiff.xi.coef function                       #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: November, 30, 2001                            #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fracdiff.xi.coef <- function(d,M) {
     xi.coef <- rep(0,M)
     secondo.termine <- (1+d*digamma(1-d))/gamma(1-d)    
     for (j in 1:M) {
          primo.termine <- gamma(j-d)/gamma(j+1)
          if (is.nan(primo.termine)) {
              primo.termine <- j^(-(1+d))*exp(d)
          }
          xi.coef[j] <- - primo.termine*(secondo.termine+digamma(j-d)/gamma(-d))
     }
return(xi.coef)
}

#############################################################
#                                                           #
#	wle.fracdiff.residuals function                     #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 11, 2001                            #
#	Version: 0.1-1                                      #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fracdiff.residuals <- function(d, M, x, x.ao, x.init=rep(0,M), x.mean=0, use.init=FALSE) {

    x <- x - x.mean
    x.ao <- x.ao - x.mean

    if (use.init) {
        MM <- 0
    } else {
        MM <- M
    }
    nused <- length(x) 
    pi.coef <- wle.fracdiff.pi.coef(d,M)
    y <- c(x.init,x.ao)
    resid <- rep(0,nused)

    for (t in (1+MM):nused) {
         resid[t] <- x[t]+pi.coef%*%y[(t-1+M):t]
    }
    resid <- resid[(MM+1):nused]

    return(resid)
}

#############################################################
#                                                           #
#	wle.fracdiff.sigma2 function                        #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 1, 2001                             #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fracdiff.sigma2 <- function(resid, w=rep(1,length(resid))) {
    sigma2 <- sum(w*resid^2)/(sum(w) - 1)
    return(sigma2)
}

#############################################################
#                                                           #
#	wle.fracdiff.fitted function                        #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 5, 2001                             #
#	Version: 0.1-1                                      #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fracdiff.fitted <- function(t, d, M, x, x.init=rep(0,M), x.mean=0, use.init=FALSE) {
 
    x <- x - x.mean

    if (use.init) {
        MM <- 0
    } else {
        MM <- M
    }
    nused <- length(x)
    pi.coef <- wle.fracdiff.pi.coef(d,M)
    y <- c(x.init,x)
    return((-pi.coef%*%y[(t-1+M):t]+x.mean))
}







