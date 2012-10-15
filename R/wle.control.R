#############################################################
#                                                           #
#	wle.glm.control function                            #
#	Author: Claudio Agostinelli and Fatemah Alqallaf    #
#	E-mail: claudio@unive.it                            #
#	Date: March, 12, 2010                               #
#	Version: 0.3                                        #
#                                                           #
#	Copyright (C) 2010 Claudio Agostinelli              #
#                          and Fatemah Alqallaf             #
#                                                           #
#############################################################

wle.glm.control <- function(boot=30, group=NULL, num.sol=1, raf=c('GKL', 'PWD', 'HD', 'NED', 'SCHI2'), tau=0.1, cutpoint=0, powerdown=1, delta=NULL, smooth=NULL, asy.smooth=0.031, tol=10^(-6), equal=10^(-3), max.iter=500, window.size=NULL, use.asymptotic=NULL, use.smooth=TRUE, mle.dispersion=FALSE, verbose=FALSE) {  
  if (!is.numeric(boot) || boot <= 0) 
    stop("value of 'boot' must be > 0")
  if(!is.null(group))
    if(!is.numeric(group) || group <= 0)
      stop("value of 'group' must be > 0 or NULL")
  if (!is.numeric(num.sol) || num.sol < 1 || !isTRUE(all.equal(as.integer(num.sol), num.sol)))
    stop("value of 'num.sol' must be an integer")
  raf <- match.arg(raf)
  if ((!is.numeric(tau) || tau > 1 || tau < 0) && raf=='GKL')
    stop("value of 'tau' must be in [0,1] when 'raf=GKL'")
  ##park+basu+2003.pdf
  if ((!is.numeric(tau) || tau < -1) && raf=='PWD')
    stop("value of 'tau' must be in [-1,Inf] when 'raf=PWD'")
  if (!is.numeric(cutpoint) || cutpoint > 1 || cutpoint < 0)
    stop("value of 'cutpoint' must be in [0,1]")
  if (!is.numeric(powerdown))
    stop("value of 'powerdown' must be numeric")  
  ##lindsay+1994.pdf
  if(!is.null(delta))
    if(!is.numeric(delta) || delta <= 0 || delta >=1)
      stop("value of 'delta' must be > 0 and < 1")
  if(!is.null(smooth))
    if(!is.numeric(smooth) || smooth <= 0)
      stop("value of 'smooth' must be > 0 or NULL")
  if(!is.numeric(asy.smooth) || asy.smooth <= 0)
    stop("value of 'asy.smooth' must be > 0")
  if (!is.numeric(tol) || tol <= 0)
    stop("value of 'tol' must be > 0")
  if (!is.numeric(equal) || equal <= tol)
    stop("value of 'equal' must be greater then 'tol'")
  if (!is.numeric(max.iter) || max.iter <= 0) 
    stop("maximum number of iterations must be > 0")
  if(!is.null(window.size))
    if(!is.numeric(window.size) || window.size <= 0)
      stop("value of 'window.size' must be > 0 or NULL")
  if(!is.null(use.asymptotic))
    if(!is.numeric(use.asymptotic) || length(use.asymptotic)!=1 || use.asymptotic <= 0)
      stop("value of 'use.asymptotic' must be a scalar > 0")
  if (!(is.logical(use.smooth) && length(use.smooth)==1))
    stop("'use.smooth' must be logical")
  if (!(is.logical(mle.dispersion) && length(mle.dispersion)==1))
    stop("'mle.dispersion' must be logical")  
  if (!(is.logical(verbose) && length(verbose)==1))
    stop("'verbose' must be logical")  
  list(boot=boot, group=group, num.sol=num.sol, raf=raf, tau=tau, cutpoint=cutpoint, powerdown=powerdown, delta=delta, smooth=smooth, asy.smooth=asy.smooth, tol=tol, equal=equal, max.iter=max.iter, window.size=window.size, use.asymptotic=use.asymptotic, use.smooth=use.smooth, mle.dispersion=mle.dispersion, verbose=verbose)
}

#############################################################
#                                                           #
#	wle.lm.control function                             #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: May, 18, 2010                                 #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2010 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.lm.control <- function(nstart=30, group=NULL, num.sol=1, raf=c('HD', 'NED', 'SCHI2', 'GKL', 'PWD'), tau=0.1, cutpoint=0, powerdown=1, smooth=0.031, tol=10^(-6), equal=10^(-3), max.iter=500, verbose=FALSE) {
  
  if (!is.numeric(nstart) || nstart <= 0) 
    stop("value of 'nstart' must be > 0")
  if(!is.null(group))
    if(!is.numeric(group) || group <= 0)
      stop("value of 'group' must be > 0 or NULL")
  if (!is.numeric(num.sol) || num.sol < 1 || !isTRUE(all.equal(as.integer(num.sol), num.sol)))
    stop("value of 'num.sol' must be an integer")
  raf <- match.arg(raf)
  if (raf=='GKL' | raf=='PWD')
    stop("no implemented for GKL or PWD RAF for now!")
  if ((!is.numeric(tau) || tau > 1 || tau < 0) && raf=='GKL')
    stop("value of 'tau' must be in [0,1] when 'raf=GKL'")
  ##park+basu+2003.pdf
  if ((!is.numeric(tau) || tau < -1) && raf=='PWD')
    stop("value of 'tau' must be in [-1,Inf] when 'raf=PWD'")
  if (!is.numeric(cutpoint) || cutpoint > 1 || cutpoint < 0)
    stop("value of 'cutpoint' must be in [0,1]")
  if (!is.numeric(powerdown))
    stop("value of 'powerdown' must be numeric")  
  ##lindsay+1994.pdf
  if(!is.numeric(smooth) || smooth <= 0)
    stop("value of 'smooth' must be > 0")
  if (!is.numeric(tol) || tol <= 0)
    stop("value of 'tol' must be > 0")
  if (!is.numeric(equal) || equal <= tol)
    stop("value of 'equal' must be greater then 'tol'")
  if (!is.numeric(max.iter) || max.iter <= 0) 
    stop("maximum number of iterations must be > 0")
  if (!(is.logical(verbose) && length(verbose)==1))
    stop("'verbose' must be logical")  
  list(nstart=nstart, group=group, num.sol=num.sol, raf=raf, tau=tau, cutpoint=cutpoint, powerdown=powerdown, smooth=smooth, tol=tol, equal=equal, max.iter=max.iter, verbose=verbose)
}
