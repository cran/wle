#############################################################
#                                                           #
#	WLE.AIC function                                    #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: October, 10, 1999                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 1999 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.aic_function(ydata,xdata,boot=100,group,inter=1,var.full=0,num.sol=1,raf=1,smooth=0.0320018,tol=10^(-6),equal=10^(-3),max.iter=500,min.weight=0.5,type=0,alpha=2)
{
size_length(ydata)
nvar_length(xdata)/length(ydata)
nrep_(2^(nvar+inter))-1

if(size<nvar+inter+2){stop("wle.aic: Number of observation must be at least equal to the number of predictors (including intercept) + 2")}
if(nvar<1){stop("Warning: wle.aic: the number of the predictors must be at least one")}
maxboot_sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))
if(boot<1 | log(boot) > maxboot){
stop("wle.aic: bootstrap replication not in the range")
}
if(!(inter==1)){inter_0}
if(!(group>1)){
group_max(round(size/4),nvar+inter+2)
warning("wle.aic: dimension of the subsample set to default value")
}
if(!(num.sol>=1)){
warning("wle.aic: number of solution to report set to 1")
num.sol_1
}
if(!(raf==1 | raf==2 | raf==3)){
warning("wle.aic: Helliger Residual Adjustment Function is used")
raf_1
}
if(max.iter<1){
warning("wle.aic: max number of iteration set to 500")
max.iter_500
}
if(smooth<10^(-5)){
warning("wle.aic: the smooth parameter seems too small")
}
if(tol<0){
warning("wle.aic: the accuracy can not be negative, using default value")
tol_10^(-6)
}
if(equal<0){
warning("wle.aic: the equal parameter can not be negative, using default value")
equal_10^(-3)
}
if(var.full<0){
warning("wle.aic: the variance of the full model can not be negative, using default value")
var.full_0
}
if(min.weight<0){
warning("wle.aic: the minimum sum of the weights can not be negative, using default value")
min.weight_0.5
}
if(!(type==0 | type==1)){
warning("wle.aic: the type must be 0 or 1, the default value is 0")
type_0
}

  z <- .Fortran("wleaic",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(inter), 
	as.integer(size),
	as.integer(nvar),
	as.integer(boot),
	as.integer(group),
	as.integer(nrep),
	as.integer(raf),
	as.double(smooth),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	as.double(var.full),
	as.integer(num.sol),
	as.double(min.weight),
	as.integer(type),
	as.double(alpha),
	cp=mat.or.vec(nrep*num.sol,nvar+inter+1),
	param=mat.or.vec(nrep*num.sol,nvar+inter),
	var=double(nrep*num.sol),
	resid=mat.or.vec(nrep*num.sol,size),
	totweight=double(nrep*num.sol),
	weight=mat.or.vec(nrep*num.sol,size),
	same=integer(nrep*num.sol),
	info=integer(1))

delnull_z$same==0
cp_z$cp[!delnull,]
param_z$param[!delnull,]
var_z$var[!delnull]
resid_z$resid[!delnull,]
totweight_z$totweight[!delnull]
weight_z$weight[!delnull,]
same_z$same[!delnull]

return(list(waic=cp,coefficients=param,scale=sqrt(var),residuals=resid,tot.weights=totweight,weights=weight,freq=same,info=z$info))

}
