#############################################################
#                                                           #
#	WLE.NORMAL.MULTI function                           #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: October, 10, 1999                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 1999 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.normal.multi_function(data,boot,group,num.sol=1,raf=1,smooth=0.00300996,tol=10^(-6),equal=10^(-3),max.iter=500)
{

data_as.matrix(data)
size_dim(data)[1]
nvar_dim(data)[2]

if(size<nvar*(nvar+1)){stop("wle.normal.multi: Number of observation must be at least equal to ")}
if(!(group>nvar*(nvar+1))){
group_max(round(size/4),nvar*(nvar+1))
warning("wle.normal.multi: dimension of the subsample set to default value")
}
maxboot_sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))
if(boot<1 | log(boot) > maxboot){
stop("wle.normal.multi: bootstrap replication not in the range")
}
if(!(num.sol>=1)){
warning("wle.normal.multi: number of solution to report set to 1")
num.sol_1
}
if(!(raf==1 | raf==2 | raf==3)){
warning("wle.normal.multi: Helliger Residual Adjustment Function is used")
raf_1
}
if(max.iter<1){
warning("wle.normal.multi: max number of iteration set to 500")
max.iter_500
}
if(smooth<10^(-5)){
warning("wle.normal.multi: the smooth parameter seems too small")
}
if(tol<0){
warning("wle.normal.multi: the accuracy can not be negative, using default value")
tol_10^(-6)
}
if(equal<0){
warning("wle.normal.multi: the equal parameter can not be negative, using default value")
equal_10^(-3)
}

  z <- .Fortran("wlenormmulti",
	as.double(data), 
	as.integer(size),
	as.integer(nvar),
	as.integer(boot),
	as.integer(group),
	as.integer(num.sol),
	as.integer(raf),
	as.double(smooth),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	mean=mat.or.vec(num.sol,nvar),
	var=mat.or.vec(num.sol,nvar*nvar),
	totweight=double(num.sol),
	weight=mat.or.vec(num.sol,size),
	same=integer(num.sol),
	nsol=integer(1),
	nconv=integer(1))



if(nvar>1)
{
temp_z$var[1:z$nsol,]
if(z$nsol>1)
{
temp.a_matrix(temp[1,],ncol=nvar)
scale_list(temp.a)
for(i in 2:z$nsol){
temp.a_matrix(temp[i,],ncol=nvar)
scale_c(scale,list(temp.a))
}
}
else
{
temp.a_matrix(temp,ncol=nvar)
scale_list(temp.a)
}
return(list(location=z$mean[1:z$nsol,],variance=scale,tot.weights=z$totweight[1:z$nsol],weights=z$weight[1:z$nsol,],freq=z$same[1:z$nsol],tot.sol=z$nsol,not.conv=z$nconv))
}
else
{
return(list(location=z$mean[1:z$nsol],variance=z$var[1:z$nsol],tot.weights=z$totweight[1:z$nsol],weights=z$weight[1:z$nsol,],freq=z$same[1:z$nsol],tot.sol=z$nsol,not.conv=z$nconv))
}

}




