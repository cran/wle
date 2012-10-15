#############################################################
#                                                           #
#	mle.cv.rb function                                  #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: May, 18, 2010                                 #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2010 Claudio Agostinelli              #
#                                                           #
#############################################################

mle.cv.rb <- function(formula, data=list(), model=TRUE, x=FALSE, y=FALSE, contrasts=NULL, monte.carlo=500, split=NULL, R=1000, sim=c('betaR', 'uniformR', 'fastRBwc'), seed=1234, num.max=NULL, nmodel=NULL, file.name=NULL, append=FALSE, mean=0.4, shape=NULL, nsize=NULL, replace=TRUE, verbose=FALSE, ...) {

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

  info <- mcv <- mcv0 <- vector(length=0)
  ret.x <- x
  ret.y <- y
  result <- list()	
  mt <- terms(formula, data = data)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$monte.carlo <- mf$split <- NULL
  mf$contrasts <- NULL
  mf$model <- mf$x <- mf$y <- NULL
  mf$verbose <- NULL
  mf$num.max <- mf$nmodel <- mf$seed <- mf$file <- mf$... <- NULL
  mf$R <- mf$sim <- mf$mean <- mf$shape <- NULL
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
    stop("mle.cv.rb: number of observations must be at least equal to the number of predictors (including intercept)")

  if (is.null(split) || split<nvar+2 || split>(size-2)) {
    split <- max(round(size^(3/4)),nvar+2)
    if (verbose)
      cat("mle.cv.rb: dimension of the split subsample set to default value = ",split,"\n")
  }

  maxcarlo <- sum(log(1:size))-(sum(log(1:split))+sum(log(1:(size-split))))

  if (monte.carlo<1 | log(monte.carlo) > maxcarlo)
    stop("mle.cv.rb: MonteCarlo replication not in the range")

  if (is.null(shape))
    shape <- length(y)
  shape1 <- shape*mean/(1-mean)  

  if (is.null(nsize))
    nsize <- size

  posizione <- 1:nvar

####################### t0
  for (imodel in nmodel) {
    if (verbose)
      cat("mle.cv.rb: t0, Running model with number ", imodel, " on ", length(nmodel), " models\n")
    zcv <- .Fortran("wlecvone",
      as.double(ydata),
      as.matrix(xdata),
      as.integer(0), 
      as.integer(size),
      as.integer(nvar),
      as.integer(monte.carlo),
      as.integer(imodel),
      as.integer(split),
      as.integer(seed),
      as.double(rep(1, size)),
      mcv=double(nvar+1),
      info=integer(1),
      PACKAGE="wle")
    info <- c(info, zcv$info)
      if (!is.null(zcv$mcv))
        mcv0 <- rbind(mcv0, zcv$mcv)
  }
  modelli <- mcv0[,1:nvar]
  t0.cv <- mcv0[,nvar+1]
  names(t0.cv) <- nmodel
  mcv0 <- as.vector(mcv0[which.min(mcv0[,nvar+1]),])
  if (NROW(mcv0)) {
    var.in <- posizione*mcv0[posizione]
    var.in <- var.in[-var.in!=0]
    t0.param <- rep(0, nvar)
    temp.lm <- summary(lm(ydata~xdata[,var.in] -1, weights=rep(1, size)))
    t0.param[var.in] <- c(temp.lm$coeff[,1])
    t0.param <- c(t0.param, temp.lm$sigma, mcv0[nvar+1])
    names(t0.param) <- c(colnames(xdata), 'scale', 'mcv')
    if (!is.null(file.param))
      write.table(x=matrix(t0.param, nrow=1), file=file.param, append=TRUE, ...)
  }

##########################  t
  if (is.null(file.name)) {
    t.param <- vector(length=0)
    t.cv <- matrix(0, nrow=R, ncol=nrep)
  }

  iboot <- 1
  while (iboot <= R) {
    if (sim=='uniformR') {
      nsize <- size
      dpesiboot <- runif(n=size, min=0, max=1)
    } else if (sim=='betaR') {
      nsize <- size
      dpesiboot <- rbeta(size, shape1=shape1, shape2=shape, ncp=0)
    } else {
      isample <- sample(x=1:size, size=nsize, replace=TRUE)
      yboot <- ydata[isample]
      xboot <- xdata[isample,]
      dpesiboot <- rep(1, nsize)
    }
    dpesiboot <- dpesiboot/sum(dpesiboot)
    mcv <- vector(length=0)
    for (imodel in nmodel) {
      if (verbose)
        cat("Running model with number ", imodel, " on ", length(nmodel), " models\n")

      zcv <- .Fortran("wlecvone",
        as.double(yboot),
        as.matrix(xboot),
        as.integer(0), 
        as.integer(nsize),
        as.integer(nvar),
        as.integer(monte.carlo),
        as.integer(imodel),
        as.integer(split),
        as.integer(seed),
        as.double(dpesiboot),
        mcv=double(nvar+1),
        info=integer(1),
        PACKAGE="wle")
      info <- c(info, zcv$info)
      if (!is.null(zcv$mcv))
        mcv <- rbind(mcv, zcv$mcv)
#  fine nmodel
    }
    if (is.null(file.cv))
      t.cv[iboot,] <- mcv[,nvar+1]
    else
      write.table(x=matrix(mcv[,nvar+1], nrow=1), file=file.cv, append=TRUE, ...)
    mcv <- as.vector(mcv[which.min(mcv[,nvar+1]),])
   
    if (NROW(mcv)) {
      iboot <- iboot + 1
      var.in <- posizione*mcv[posizione]
      var.in <- var.in[-var.in!=0]
      param <- rep(0, nvar)
      temp.lm <- summary(lm(ydata~xdata[,var.in] -1,weights=dpesiboot))
      param[var.in] <- temp.lm$coeff[,1]
      param <- c(param, temp.lm$sigma, mcv[nvar+1])

      if (is.null(file.param))
        t.param <- rbind(t.param, param)
      else
        write.table(x=matrix(param, nrow=1), file=file.param, append=TRUE, ...)
    }
# fine bootstrap
  }

  colnames(t.param) <- c(colnames(xdata), 'scale', 'mcv')
  colnames(modelli) <- colnames(xdata)
  rownames(modelli) <- colnames(t.cv) <- nmodel
  rownames(t.param) <- rownames(t.cv) <- 1:R
  
  if (is.null(file.name)) {
    result$t0 <- t0.param
    result$t0.cv <- t0.cv
    result$t <- t.param
    result$t.cv <- t.cv    
    result$models.matrix <- modelli
    result$call <- cl
    result$info <- info
    result$contrasts <- attr(xdata, "contrasts")
    result$xlevels <- xlev
    result$terms <- mt
    result$seed <- seed
    if (model)
        result$model <- mf
    if (ret.x)
        result$x <- xdata
    if (ret.y)
        result$y <- ydata
    class(result) <- 'mle.cv.rb'
  } else {
    result <- paste("Results saved to files ", file.param, " and to ", file.cv, "\n", sep='')
    if (verbose)
      cat(result)
  }
  return(result)
}

