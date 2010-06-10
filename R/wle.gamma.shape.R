#############################################################
#                                                           #
#	wle.glm function                                    #
#	Author: Claudio Agostinelli and Fatemah Alqallaf    #
#	E-mail: claudio@unive.it                            #
#	Date: December, 15, 2009                            #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2009 Claudio Agostinelli              #
#                      and Fatemah Alqallaf                 #
#                                                           #
#############################################################

## Function developed from 'MASS/R/gamma.shape.R' version 7.2-48
# copyright (C) 1994-2009 W. N. Venables and B. D. Ripley

wle.gamma.shape.glm <- function(y, mu, deviance, df.residual, prior.weights=NULL, wle.weights=NULL, it.lim = 10, eps.max = .Machine$double.eps^0.25, verbose = FALSE, ...) {
    if (is.null(prior.weights))
      prior.weights <- rep(1, length(y))
    if (any(prior.weights!=1))
      stop("In the 'Gamma' family 'prior.weigths' must be equal to ones, no implementation is available for the general case")
    if (is.null(wle.weights))
      wle.weights <- rep(1, length(y))    
    Dbar <- deviance/df.residual
    alpha <- (6 + 2*Dbar)/(Dbar*(6 + Dbar))
    if(verbose) {
	message("Initial estimate: ", format(alpha))
	utils::flush.console()
    }
    fixed <-  -y/mu - log(mu) + log(prior.weights) + 1 + log(y + (y == 0))
    eps <- 1
    itr <- 0
    while(abs(eps) > eps.max && (itr <- itr + 1) <= it.lim) {
        sc <- sum(wle.weights * prior.weights * (fixed + log(alpha) - digamma(prior.weights * alpha)))
        inf <- sum(wle.weights * prior.weights * (prior.weights * trigamma(prior.weights * alpha) - 1/alpha))
        alpha <- alpha + (eps <- sc/inf)
        if(verbose) {
	    message("Iter. ", itr, " Alpha: ", format(alpha))
	    utils::flush.console()
	}
    }
    if(itr > it.lim) warning("iteration limit reached")
    res <- list(dispersion=1/alpha, alpha = alpha, SE = sqrt(1/inf))
    return(res)
}


