#############################################################
#                                                           #
#	MLE.CV function                                     #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: October, 10, 1999                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 1999 Claudio Agostinelli              #
#                                                           #
#############################################################

mle.cv_function(ydata,xdata,monte.carlo=500,split,inter=1)
{
size_length(ydata)
nvar_length(xdata)/length(ydata)
nrep_(2^(nvar+inter))-1

if(size<nvar+inter+2){stop("mle.cv: Number of observation must be at least equal to the number of predictors (including intercept) + 2")}
if(nvar<1){stop("mle.cv: the number of the predictors must be at least one")}

if(!(inter==1)){inter_0}

if(split<nvar+inter+2 | split>(size-2)){
split_max(round(size^(3/4)),nvar+inter+2)
warning(paste("mle.cv: dimension of the split subsample set to default value = ",split))
}
maxcarlo_sum(log(1:size))-(sum(log(1:split))+sum(log(1:(size-split))))
if(monte.carlo<1 | log(monte.carlo) > maxcarlo){
stop("mle.cv: MonteCarlo replication not in the range")
}

  z <- .Fortran("mlecv",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(inter), 
	as.integer(size),
	as.integer(nvar),
	as.integer(nrep),
	as.integer(monte.carlo),
	as.integer(split),
	cv=mat.or.vec(nrep,nvar+inter+1),
	info=integer(1))

return(list(cv=z$cv,info=z$info))

}

