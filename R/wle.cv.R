#############################################################
#                                                           #
#	WLE.CV function                                     #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: October, 10, 1999                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 1999 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.cv_function(ydata,xdata,boot=100,group,monte.carlo=500,split,inter=1,num.sol=1,raf=1,smooth=0.0320018,tol=10^(-6),equal=10^(-3),max.iter=500,min.weight=0.5)
{
size_length(ydata)
nvar_length(xdata)/length(ydata)
nrep_(2^(nvar+inter))-1

if(size<nvar+inter+2){stop("wle.cv: Number of observation must be at least equal to the number of predictors (including intercept) + 2")}
if(nvar<1){stop("wle.cv: the number of the predictors must be at least one")}

if((!(group>1))|group<(nvar+inter+1)){
group_max(round(size/4),nvar+inter+2)
warning("wle.lm: dimension of the subsample set to default value")
}

maxboot_sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))
if(boot<1 | log(boot) > maxboot){
stop("wle.cv: bootstrap replication not in the range")
}
if(!(inter==1)){inter_0}

if(split<nvar+inter+2 | split>(size-2)){
split_max(round(size^(3/4)),nvar+inter+2)
warning(paste("wle.cv: dimension of the split subsample set to default value = ",split))
}
maxcarlo_sum(log(1:size))-(sum(log(1:split))+sum(log(1:(size-split))))
if(monte.carlo<1 | log(monte.carlo) > maxcarlo){
stop("wle.cv: MonteCarlo replication not in the range")
}
if(!(num.sol>=1)){
warning("wle.cv: number of solution to report set to 1")
num.sol_1
}
if(!(raf==1 | raf==2 | raf==3)){
warning("wle.cv: Helliger Residual Adjustment Function is used")
raf_1
}
if(max.iter<1){
warning("wle.cv: max number of iteration set to 500")
max.iter_500
}
if(smooth<10^(-5)){
warning("wle.cv: the smooth parameter seems too small")
}
if(tol<0){
warning("wle.cv: the accuracy can not be negative, using default value")
tol_10^(-6)
}
if(equal<0){
warning("wle.cv: the equal parameter can not be negative, using default value")
equal_10^(-3)
}
if(min.weight<0){
warning("wle.cv: the minimum sum of the weights can not be negative, using default value")
min.weight_0.5
}



  z <- .Fortran("wlecv",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(inter), 
	as.integer(size),
	as.integer(nvar),
	as.integer(boot),
	as.integer(group),
	as.integer(nrep),
	as.integer(monte.carlo),
	as.integer(split),
	as.integer(raf),
	as.double(smooth),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	as.integer(num.sol),
	as.double(min.weight),
	cv=mat.or.vec(nrep,nvar+inter+1),
	param=mat.or.vec(num.sol,nvar+inter),
	var=double(num.sol),
	resid=mat.or.vec(num.sol,size),
	totweight=double(num.sol),
	weight=mat.or.vec(num.sol,size),
	same=integer(num.sol),
	index=integer(1),
	info=integer(1))

delnull_z$same==0
param_z$param[!delnull,]
var_z$var[!delnull]
resid_z$resid[!delnull,]
totweight_z$totweight[!delnull]
weight_z$weight[!delnull,]
same_z$same[!delnull]

return(list(wcv=z$cv,coefficients=param,scale=sqrt(var),residuals=resid,tot.weights=totweight,weights=weight,freq=same,index=z$index,info=z$info))

}

