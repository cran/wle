#############################################################
#                                                           #
#	WLE.STEPWISE function                               #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: October, 10, 1999                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 1999 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.stepwise_function(ydata,xdata,boot=100,group,inter=1,num.sol=1,raf=1,smooth=0.0320018,tol=10^(-6),equal=10^(-3),max.iter=500,min.weight=0.5,type="Forward",f.in=0.0,f.out=0.0,method="WLE")
{
size_length(ydata)
xdata_as.matrix(xdata)
nvar_length(c(xdata))/length(ydata)
nrep_(2^(nvar+inter))-1

if(size<nvar+inter+2){stop("wle.stepwise: Number of observation must be at least equal to the number of predictors (including intercept) + 2")}
if(nvar<1){stop("wle.stepwise: The number of the predictors must be at least one")}
maxboot_sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))
if(boot<1 | log(boot) > maxboot){
stop("wle.stepwise: bootstrap replication not in the range")
}
if(!(inter==1)){inter_0}
if(!(group>1)){
group_max(round(size/4),nvar+inter+2)
warning("wle.stepwise: dimension of the subsample set to default value")
}
if(!(num.sol>=1)){
warning("wle.stepwise: number of solution to report set to 1")
num.sol_1
}
if(!(raf==1 | raf==2 | raf==3)){
warning("wle.stepwise: Helliger Residual Adjustment Function is used")
raf_1
}
if(max.iter<1){
warning("wle.stepwise: max number of iteration set to 500")
max.iter_500
}
if(smooth<10^(-5)){
warning("wle.stepwise: the smooth parameter seems too small")
}
if(tol<0){
warning("wle.stepwise: the accuracy can not be negative, using default value")
tol_10^(-6)
}
if(equal<0){
warning("wle.stepwise: the equal parameter can not be negative, using default value")
equal_10^(-3)
}
if(min.weight<0){
warning("wle.stepwise: the minimum sum of the weights can not be negative, using default value")
min.weight_0.5
}
if(!(type=="Forward" | type=="Backward" | type=="Stepwise")){
warning("wle.stepwise: the type must be Forward, Backward or Stepwise, the default value is Forward")
type_"Forward"
}

if(type=="Forward"){ntype_1}
if(type=="Backward"){ntype_2}
if(type=="Stepwise"){ntype_3}

if(!(method=="WLE" | method=="WLS")){
warning("wle.stepwise: the method must be WLE, or WLS the default value is WLE")
method_"WLE"

}

if(method=="WLE"){nmethod_0}
if(method=="WLS"){nmethod_1}

  z <- .Fortran("wstep",
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
	as.integer(ntype),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	as.integer(num.sol),
	as.double(min.weight),
	as.double(f.in),
	as.double(f.out),
	as.integer(nmethod),
	step=mat.or.vec(nrep,nvar+inter+1),
	param=mat.or.vec(num.sol,nvar+inter),
	var=double(num.sol),
	resid=mat.or.vec(num.sol,size),
	totweight=double(num.sol),
	weight=mat.or.vec(num.sol,size),
	same=integer(num.sol),
	indice=integer(1),
	info=integer(1),
	imodel=integer(1),
	nsol=integer(1))

step_z$step[1:z$imodel,]
param_z$param[1:z$nsol,]
var_z$var[1:z$nsol]
resid_z$resid[1:z$nsol,]
totweight_z$totweight[1:z$nsol]
weight_z$weight[1:z$nsol,]
same_z$same[1:z$nsol]

return(list(wstep=step,coefficients=param,scale=sqrt(var),residuals=resid,tot.weights=totweight,weights=weight,tot.sol=z$nsol,freq=same,index=z$indice,info=z$info))

}
