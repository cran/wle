#############################################################
#                                                           #
#	WLE.NORMAL function                                 #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: October, 10, 1999                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 1999 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.normal_function(data,boot,group,num.sol=1,raf=1,smooth=0.00300996,tol=10^(-6),equal=10^(-3),max.iter=500)
{

data_as.vector(data)
size_length(data)

if(size<3){stop("wle.normal: Number of observation must be at least equal to 3")}
if(!(group>1)){
group_max(round(size/4),3)
warning("wle.normal: dimension of the subsample set to default value")
}
maxboot_sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))
if(boot<1 | log(boot) > maxboot){
stop("wle.normal: bootstrap replication not in the range")
}
if(!(num.sol>=1)){
warning("wle.normal: number of solution to report set to 1")
num.sol_1
}
if(!(raf==1 | raf==2 | raf==3)){
warning("wle.normal: Helliger Residual Adjustment Function is used")
raf_1
}
if(max.iter<1){
warning("wle.normal: max number of iteration set to 500")
max.iter_500
}
if(smooth<10^(-5)){
warning("wle.normal: the smooth parameter seems too small")
}
if(tol<0){
warning("wle.normal: the accuracy can not be negative, using default value")
tol_10^(-6)
}
if(equal<0){
warning("wle.normal: the equal parameter can not be negative, using default value")
equal_10^(-3)
}

  z <- .Fortran("wlenorm",
	as.double(data), 
	as.integer(size),
	as.integer(boot),
	as.integer(group),
	as.integer(num.sol),
	as.integer(raf),
	as.double(smooth),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	mean=double(num.sol),
	var=double(num.sol),
	totweight=double(num.sol),
	weight=mat.or.vec(num.sol,size),
	same=integer(num.sol),
	nsol=integer(1),
	nconv=integer(1))





return(list(location=z$mean[1:z$nsol],scale=sqrt(z$var)[1:z$nsol],tot.weights=z$totweight[1:z$nsol],weights=z$weight[1:z$nsol,],freq=z$same[1:z$nsol],tot.sol=z$nsol,not.conv=z$nconv))

}




