#############################################################
#                                                           #
#	WLE.LM function                                     #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: October, 10, 1999                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 1999 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.lm_function(ydata,xdata,boot=100,group,inter=1,num.sol=1,raf=1,smooth=0.0320018,tol=10^(-6),equal=10^(-3),max.iter=500)
{

xdata_as.matrix(xdata)

size_length(ydata)
nvar_length(xdata)/length(ydata)
ncol_nvar

if(size<(nvar+inter+2)){stop("wle.lm: Number of observation must be at least equal to the number of predictors (including intercept) + 2")}
if(nvar<1){stop("wle.lm: the number of the predictors must be at least one")}

if((!(group>1))|group<(nvar+inter+1)){
group_max(round(size/4),nvar+inter+2)
warning("wle.lm: dimension of the subsample set to default value")
}

maxboot_sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))
if(boot<1 | log(boot) > maxboot){
stop("wle.lm: bootstrap replication not in the range")
}
if(!(inter==1)){inter_0}
if(!(num.sol>=1)){
warning("wle.lm: number of solution to report set to 1")
num.sol_1
}
if(!(raf==1 | raf==2 | raf==3)){
warning("wle.lm: Helliger Residual Adjustment Function is used")
raf_1
}
if(max.iter<1){
warning("wle.lm: max number of iteration set to 500")
max.iter_500
}
if(smooth<10^(-5)){
warning("wle.lm: the smooth parameter seems too small")
}
if(tol<0){
warning("wle.lm: the accuracy can not be negative, using default value")
tol_10^(-6)
}
if(equal<0){
warning("wle.lm: the equal parameter can not be negative, using default value")
equal_10^(-3)
}

  z <- .Fortran("wleregfix",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(inter), 
	as.integer(size),
	as.integer(ncol),
	as.integer(nvar),
	as.integer(boot),
	as.integer(group),
	as.integer(num.sol),
	as.integer(raf),
	as.double(smooth),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	param=mat.or.vec(num.sol,nvar+inter),
	var=double(num.sol),
	resid=mat.or.vec(num.sol,size),
	totweight=double(num.sol),
	weight=mat.or.vec(num.sol,size),
	same=integer(num.sol),
	nsol=integer(1),
	nconv=integer(1))

if(inter==1){ xidata_cbind(xdata,rep(1,size)) } else { xidata_xdata }

if(z$nsol>0)
{
z$var_z$var[1:z$nsol]
z$totweight_z$totweight[1:z$nsol]
z$same_z$same[1:z$nsol]

if(num.sol==1)
{
z$param_c(z$param)
z$resid_c(z$resid)
z$weight_c(z$weight)
}
else
{
z$param_z$param[1:z$nsol,]
z$resid_z$resid[1:z$nsol,]
z$weight_z$weight[1:z$nsol,]
}

y.fit_t(xidata%*%matrix(z$param,ncol=z$nsol,byrow=T))

if(z$nsol==1)
{
devparam_sqrt(z$var*diag(solve(t(xidata)%*%diag(z$weight)%*%xidata,tol=1e-100)))
y.fit_c(y.fit)
}
else
{
devparam_sqrt(z$var[1]*diag(solve(t(xidata)%*%diag(z$weight[1,])%*%xidata,tol=1e-100)))
for(i in 2:z$nsol){
devparam_rbind(devparam,sqrt(z$var[i]*diag(solve(t(xidata)%*%diag(z$weight[i,])%*%xidata,tol=1e-100))))
}
}


return(list(coefficients=z$param,standard.error=devparam,scale=sqrt(z$var),tot.weights=z$totweight,freq=z$same,tot.sol=z$nsol,not.conv=z$nconv,residuals=z$resid,weights=z$weight,fitted.values=y.fit,ydata=ydata,xdata=xidata))

}
else
{
stop("wle.lm: No solutions are fuond, checks the parameters")
#,not.conv=z$nconv)
}

}

