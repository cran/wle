#############################################################
#                                                           #
#	MLE.CP function                                     #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: October, 10, 1999                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 1999 Claudio Agostinelli              #
#                                                           #
#############################################################

mle.cp_function(ydata,xdata,inter=1,var.full=0)
{
size_length(ydata)
nvar_dim(xdata)[2]
nrep_(2^(nvar+inter))-1

if(size<nvar+inter+2){stop("mle.cp: Number of observation must be at least equal to the number of predictors (including intercept) + 2")}
if(nvar<1){stop("mle.cp: the number of the predictors must be at least one")}
if(!(inter==1)){inter_0}
if(var.full<0){
warning("mle.cp: the variance of the full model can not be negative, using default value")
var.full_0
}

  z <- .Fortran("mlecp",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(inter), 
	as.integer(size),
	as.integer(nvar),
	as.integer(nrep),
	as.double(var.full),
	cp=mat.or.vec(nrep,nvar+inter+1),
	param=mat.or.vec(nrep,nvar+inter),
	var=double(nrep),
	resid=mat.or.vec(nrep,size),
	info=integer(1))

	
cp_z$cp
param_z$param
var_z$var
resid_z$resid

return(list(cp=cp,coefficients=param,scale=sqrt(var),residuals=resid,info=z$info))

}


