#############################################################
#                                                           #
#	MLE.AIC function                                    #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: October, 10, 1999                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 1999 Claudio Agostinelli              #
#                                                           #
#############################################################





mle.aic_function(ydata,xdata,inter=1,var.full=0,alpha=2)
{
size_length(ydata)
nvar_dim(xdata)[2]
nrep_(2^(nvar+inter))-1

if(size<nvar+inter+2){stop("mle.aic: Number of observation must be at least equal to the number of predictors (including intercept) + 2")}
if(nvar<1){stop("mle.aic: the number of the predictors must be at least one")}
if(!(inter==1)){inter_0}
if(var.full<0){
warning("mle.aic: the variance of the full model can not be negative, using default value")
var.full_0
}

  z <- .Fortran("mleaic",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(inter), 
	as.integer(size),
	as.integer(nvar),
	as.integer(nrep),
	as.double(var.full),
	as.double(alpha),
	aic=mat.or.vec(nrep,nvar+inter+1),
	param=mat.or.vec(nrep,nvar+inter),
	var=double(nrep),
	resid=mat.or.vec(nrep,size),
	info=integer(1))

	
aic_z$aic
param_z$param
var_z$var
resid_z$resid

return(list(aic=aic,coefficients=param,scale=sqrt(var),residuals=resid,info=z$info))

}


