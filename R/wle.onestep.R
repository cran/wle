#############################################################
#                                                           #
#	WLE.ONESTEP function                                #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: October, 10, 1999                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 1999 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.onestep_function(ydata,xdata,ini.param,ini.scale,inter=1,raf=1,smooth=0.0320018,num.step=1)
{

xdata_as.matrix(xdata)

size_length(ydata)
nvar_dim(xdata)[2]
ncol_nvar

if(size<(nvar+inter+2)){stop("wle.onestep: Number of observation must be at least equal to the number of predictors (including intercept) + 2")}
if(nvar<1){stop("wle.onestep: the number of the predictors must be at least one")}
if(!(ini.scale>=0)){
stop("wle.onestep: the initial scale error must be non negative")
}
if(!(inter==1)){inter_0}
if(!(num.step>=1)){
warning("wle.onestep: number of steps can not be negative, set to 1")
num.step_1
}
if(!(raf==1 | raf==2 | raf==3)){
warning("wle.onestep: Helliger Residual Adjustment Function is used")
raf_1
}
if(smooth<10^(-5)){
warning("wle.onestep: the smooth parameter seems too small")
}

ini.var_ini.scale^2

  z <- .Fortran("wleonestepfix",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(inter), 
	as.integer(size),
	as.integer(ncol),
	as.integer(nvar),
	as.double(ini.param),
	as.double(ini.var),
	as.integer(raf),
	as.double(smooth),
	as.integer(num.step),
	param=double(nvar+inter),
	var=double(1),
	resid=double(size),
	totweight=double(1),
	weight=double(size))

if(inter==1){
xidata_cbind(xdata,rep(1,size))
}
else
{
xidata_xdata
}

devparam_sqrt(z$var*diag(solve(t(xidata)%*%diag(z$weight)%*%xidata)))

return(list(coefficients=z$param,standard.error=devparam,scale=sqrt(z$var),residuals=z$resid,tot.weights=z$totweight,weights=z$weight))
}


