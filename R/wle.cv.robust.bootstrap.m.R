#############################################################
#BOOTSTRAPPO LA MEDIA DELLE STIME DEI PARAMETRI DEL MONTECARLO #
#	wle.cv.rb function                                  #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: October, 13, 2012                             #
#	Version: 0.5                                        #
#                                                           #
#	Copyright (C) 2012 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.cv.rb.m <- function(formula, data=list(), model=TRUE, x=FALSE, y=FALSE, contrasts=NULL, monte.carlo=500, split=NULL, R=1000, sim=c('betaR', 'uniformR', 'fastRBwc'), weights=NULL, min.weight=0.5, num.max=NULL, nmodel=NULL, file.name=NULL, append=FALSE,  mean=0.4, shape=NULL, nsize=NULL, boot.resid=FALSE, boot.param=FALSE, replace=TRUE, control=wle.lm.control(), ...) {
  
  sim <- match.arg(sim)
  if (is.null(file.name)) {
    file.param <- file.cv <- NULL
  } else {
    file.param <- paste('param-', file, sep='') 
    file.cv <- paste('cv-', file, sep='')
    if (!append) {
      if (file.exists(file.param)) file.remove(file.param)
      if (file.exists(file.cv)) file.remove(file.cv)
    }
  }

  control$raf <- switch(control$raf,
    HD = 1,
    NED = 2,
    SCHI2 = 3,
    -1)

  if (control$raf==-1) stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

  info <- vector(length=0)
  ret.x <- x
  ret.y <- y
  result <- list()	
  mt <- terms(formula, data = data)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$monte.carlo <- mf$split <- NULL
  mf$weights <- NULL
  mf$min.weight <- NULL
  mf$contrasts <- NULL
  mf$model <- mf$x <- mf$y <- NULL
  mf$num.max <- mf$nmodel <- mf$file.name <- mf$... <- NULL
  mf$R <- mf$sim <- mf$mean <- mf$shape <- NULL
  mf$control <- NULL
  mf$boot.resid <- mf$boot.param <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent()))
  xvars <- as.character(attr(mt, "variables"))[-1]
  inter <- attr(mt, "intercept")
  if((yvar <- attr(mt, "response")) > 0)
    xvars <- xvars[-yvar]
  xlev <- if(length(xvars) > 0) {
    xlev <- lapply(mf[xvars], levels)
    xlev[!sapply(xlev, is.null)]
  }
  yboot <- ydata <- model.response(mf, "numeric")
  if (is.empty.model(mt)) 
    stop("The model is empty")
  else 
    xboot <- xdata <- model.matrix(mt, mf, contrasts)

  if (is.null(size <- nrow(xdata)) | is.null(nvar <- ncol(xdata)))
    stop("'x' must be a matrix")
  if (length(ydata)!=size)
    stop("'y' and 'x' are not compatible")

  nrep <- 2^nvar-1
  if (is.null(nmodel))
    nmodel <- 1:nrep
  if (is.null(num.max))
    num.max <- length(nmodel)

  if (size<nvar)
    stop("Number of observations must be at least equal to the number of predictors (including intercept)")

  if (is.null(control$group) || control$group<nvar) {
    control$group <- max(round(size/4),(nvar+1))
    control$group <- min(size, control$group)
    if (control$verbose) cat("wle.cv.rb: dimension of the subsample set to default value = ",control$group,"\n")
  }

  maxnstart <- sum(log(1:size))-(sum(log(1:control$group))+sum(log(1:(size-control$group))))

  if (control$nstart < 1 | log(control$nstart) > maxnstart)
    stop("wle.cv.rb: Bootstrap replication not in the range")

  if (is.null(split) || split<nvar+2 || split>(size-2)) {
    split <- max(round(size^(3/4)),nvar+2)
    if (control$verbose)
      cat("wle.cv.rb: dimension of the split subsample set to default value = ",split,"\n")
  }

  maxcarlo <- sum(log(1:size))-(sum(log(1:split))+sum(log(1:(size-split))))

  if (monte.carlo<1 | log(monte.carlo) > maxcarlo)
    stop("wle.cv.rb: MonteCarlo replication not in the range")
  if (!(control$num.sol>=1)) {
    if (control$verbose)
      cat("wle.cv.rb: number of solution to report set to 1 \n")
    control$num.sol <- 1
  }

  if (control$max.iter<1) {
    if (control$verbose)
      cat("wle.cv.rb: max number of iteration set to 500 \n")
    control$max.iter <- 500
  }

  if (control$smooth<10^(-5)) {
    if (control$verbose)
      cat("wle.cv.rb: the smooth parameter seems too small \n")
  }

  if (control$tol<=0) {
    if (control$verbose)
      cat("wle.cv.rb: the accuracy must be positive, using default value: 10^(-6) \n")
    control$tol <- 10^(-6)
  }

  if (control$equal<=control$tol) {
    if (control$verbose)
      cat("wle.cv.rb: the equal parameter must be greater than tol, using default value: tol+10^(-3) \n")
    control$equal <- control$tol+10^(-3)
  }

  if (min.weight<0) {
    if (control$verbose)
      cat("wle.cv: the minimum sum of the weights can not be negative, using default value \n")
    min.weight <- 0.5
  }

  if (is.null(shape))
    shape <- length(y)
  shape1 <- shape*mean/(1-mean)
  
  if (is.null(nsize))
    nsize <- size
  
  if (is.null(weights)) {
    
    z <- .Fortran("wleregfix",
      as.double(ydata),
      as.matrix(xdata),
      as.integer(0), 
      as.integer(size),
      as.integer(nvar),
      as.integer(nvar),
      as.integer(control$nstart),
      as.integer(control$group),
      as.integer(control$num.sol),
      as.integer(control$raf),
      as.double(control$smooth),
      as.double(control$tol),
      as.double(control$equal),
      as.integer(control$max.iter),
      param=mat.or.vec(control$num.sol,nvar),
      var=double(control$num.sol),
      resid=mat.or.vec(control$num.sol,size),
      totweight=double(control$num.sol),
      weight=mat.or.vec(control$num.sol,size),
      density=mat.or.vec(control$num.sol,size),
      model=mat.or.vec(control$num.sol,size),
      delta=mat.or.vec(control$num.sol,size),
      same=integer(control$num.sol),
      nsol=integer(1),
      nconv=integer(1),
      PACKAGE = "wle")

    tot.sol <- z$nsol
    if (z$nconv==control$nstart)
      stop("wle.cv.rb: No solutions in the full model")
##          info=1

    index <- 0
    wvaria <- z$var
    wtotpesi <- z$totweight
    dvar <- wvaria[1]+1
    for (i in 1:tot.sol) {
      if (dvar>wvaria[i] & min.weight<wtotpesi[i]) {
        dvar <- wvaria[i]
        index <- i
      }
    }

    if (index==0)
      stop("wle.cv.rb: No solutions in the full model satisfied the criterion")
###       info=2
    dpesi <- z$weight[index,]
    freq <- z$same[index]

    if (control$verbose)
      cat("wle.cv.rb: Found ", tot.sol, "solution/s \n")
    if (boot.resid) {
      resid <- z$resid[index,]
      coeff <- z$coefficients[index,]
      sigmafull <- z$scale[index,]
    }
  } else {
    dpesi <- weights
    index <- 1
    tot.sol <- 1
    freq <- NA
    if (boot.resid) {
      z <- lm(ydata~xdata -1, weights=dpesi)
      resid <- z$residuals
      coeff <- z$coefficients
      sigmafull <- summary(z)$sigma
    }
  }

  posizione <- 1:nvar

  wcv0 <- param0 <- scale0 <- vector(length=0)
####################### t0
  for (imodel in nmodel) {
    if (control$verbose)
      cat("wle.cv.rb: t0, Running model with number ", imodel, " on ", length(nmodel), " models\n")
    nzero <- 0
    zcv <- .Fortran("wlecvonem",
      as.double(ydata),
      as.matrix(xdata),
      as.integer(nzero), 
      as.integer(size),
      as.integer(nvar),
      as.integer(monte.carlo),
      as.integer(imodel),
      as.integer(split),
      as.double(dpesi),
      wcv=double(nvar+nzero+1),
      param=double(nvar+nzero),
      sigma=double(1),
      info=integer(1),
      PACKAGE="wle")
    info <- c(info, zcv$info)
      if (!is.null(zcv$wcv)) {
        wcv0 <- rbind(wcv0, zcv$wcv)
        param0 <- rbind(param0, zcv$param)
        scale0 <- c(scale0, zcv$sigma)
      }
  }
  modelli <- wcv0[,1:nvar]
  t0.cv <- wcv0[,nvar+1]
  names(t0.cv) <- nmodel
  if (!is.null(file.cv))
    write.table(x=matrix(t0.cv, nrow=1), file=file.cv, append=TRUE, ...)
  b0pos <- which.min(wcv0[,nvar+1])
  wcv0 <- as.vector(wcv0[b0pos,])
  param0 <- as.vector(param0[b0pos,])
  scale0 <- scale0[b0pos]
  t0.param <- c(param0, scale0, wcv0[nvar+1])
  names(t0.param) <- c(colnames(xdata), 'scale', 'wcv')
  if (!is.null(file.param))
    write.table(x=matrix(t0.param, nrow=1), file=file.param, append=TRUE, ...)

##########################  t
  if (is.null(file.name)) {
    t.cv <- matrix(0, nrow=R, ncol=nrep)
    t.param <- matrix(0, nrow=R, ncol=nvar+2)
  }
  
  iboot <- 1
  while (iboot <= R) {
    if (sim=='uniformR') {
      if (boot.resid) {
        if (boot.param) {
          yboot <- c(xdata%*%coeff+rnorm(length(ydata), mean=0, sd=sigmafull))
          dpesiboot <- runif(n=length(ydata), min=0, max=1)
        } else {
          isample <- sample(x=1:size, size=nsize, replace=replace)
          yboot <- c(xdata%*%coeff+resid[isample])
          dpesiboot <- sapply(dpesi[isample], function(x) runif(n=1, min=0, max=x))
        }  
      } else {
        dpesiboot <- sapply(dpesi, function(x) runif(n=1, min=0, max=x))
      }
    } else if (sim=='betaR') {
      if (boot.resid) {
        if (boot.param) {
          yboot <- c(xdata%*%coeff+rnorm(length(ydata), mean=0, sd=sigmafull))
          dpesiboot <- rbeta(length(dpesi), shape1=shape1, shape2=shape, ncp=0) 
        } else {
          isample <- sample(x=1:size, size=nsize, replace=replace)
          yboot <- c(xdata%*%coeff+resid[isample])
          dpesiboot <- dpesi[isample]*rbeta(length(dpesi), shape1=shape1, shape2=shape, ncp=0)
        }
      } else {
        dpesiboot <- dpesi*rbeta(length(dpesi), shape1=shape1, shape2=shape, ncp=0)
      }
    } else {
      if (boot.param) {
          yboot <- c(xdata%*%coeff+rnorm(length(ydata), mean=0, sd=sigmafull))
          dpesiboot <- rep(1, length(ydata))
      } else {
        isample <- sample(x=1:size, size=nsize, replace=replace)
        if (boot.resid) {
          yboot <- c(xdata%*%coeff+resid[isample])
        } else {
          yboot <- ydata[isample]
          xboot <- xdata[isample,]
        }
        dpesiboot <- as.vector(dpesi[isample])
      }
    }
    dpesiboot <- dpesiboot/sum(dpesiboot)*sum(dpesi)
    wcv <- param <- scale <- vector(length=0)
    for (imodel in nmodel) {
      if (control$verbose)
        cat("t: Running model with number ", imodel, " on ", length(nmodel), " models\n")
      nzero <- 0
      zcv <- .Fortran("wlecvonem",
	as.double(yboot),
	as.matrix(xboot),
	as.integer(nzero), 
	as.integer(nsize),
	as.integer(nvar),
	as.integer(monte.carlo),
        as.integer(imodel),
	as.integer(split),
        as.double(dpesiboot),
	wcv=double(nvar+nzero+1),
        param=double(nvar+nzero),
        scale=double(1),
	info=integer(1),
	PACKAGE="wle")
      info <- c(info, zcv$info)
      if (!is.null(zcv$wcv)) {
        wcv <- rbind(wcv, zcv$wcv)
        param <- rbind(param, zcv$param)
        scale <- c(scale, zcv$scale)
      }
#  fine nmodel
    }
    if (is.null(file.cv))
      t.cv[iboot,] <- wcv[,nvar+1]
    else
      write.table(x=matrix(wcv[,nvar+1], nrow=1), file=file.cv, append=TRUE, ...)

    bpos <- which.min(wcv[,nvar+1])
    wcv <- as.vector(wcv[bpos,])
    param <- param[bpos,]
    scale <- scale[bpos]
    param <- c(param, scale, wcv[nvar+1])

    if (is.null(file.param))
      t.param[iboot,] <- param
    else
      write.table(x=matrix(param, nrow=1), file=file.param, append=TRUE, ...)
    iboot <- iboot + 1 ### new bootstrap replication index
# fine bootstrap
  }

  colnames(t.param) <- c(colnames(xdata), 'scale', 'wcv')
  colnames(modelli) <- colnames(xdata)
  rownames(modelli) <- colnames(t.cv) <- nmodel
  rownames(t.param) <- rownames(t.cv) <- 1:R
  
  if (is.null(file.name)) {
    result$t0 <- t0.param
    result$t0.cv <- t0.cv
    result$t <- t.param
    result$t.cv <- t.cv
    result$weights <- dpesi
    result$tot.weights <- mean(dpesi)
    result$freq <- freq
    result$models.matrix <- modelli
    result$call <- cl
    result$info <- info
    result$index <- index
    result$tot.sol <- tot.sol
    result$contrasts <- attr(xdata, "contrasts")
    result$xlevels <- xlev
    result$terms <- mt
    if (model)
        result$model <- mf
    if (ret.x)
        result$x <- xdata
    if (ret.y)
        result$y <- ydata
    class(result) <- 'wle.cv.rb'
  } else {
    result <- paste("Results saved to files ", file.param, " and to ", file.cv, "\n", sep='')
    if (control$verbose)
      cat(result)
  }
  return(result)
}
