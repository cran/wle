#############################################################
#                                                           #
#	PLOT.WLE.LM function                                #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: October, 10, 1999                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 1999 Claudio Agostinelli              #
#                                                           #
#############################################################

plot.wle.lm_function(object,level.weight=0.5,plot.it=TRUE)
{
param_as.matrix(object$coefficients)
res_as.matrix(object$residuals)
y.fit_as.matrix(object$fitted.values)
ydata_object$ydata
xdata_object$xdata
weight_as.matrix(object$weights)
tot.weight_object$tot.weights
tot.sol_object$tot.sol

if(tot.sol>1){

x11()
par(mfcol=c(tot.sol,tot.sol))

for(isol in 1:tot.sol)
{
for(jsol in 1:tot.sol)
{
y.fit.i_y.fit[isol,]
res.i_res[isol,]
weight.i_weight[isol,]

y.fit.j_y.fit[jsol,]
res.j_res[jsol,]
weight.j_weight[jsol,]

level.i_weight.i>=level.weight
level.j_weight.j>=level.weight

color_rep(2,length(ydata))
color[level.i]_1
color[level.j]_3

color.res_rep(2,length(ydata))
color.res[level.i & (res.i>res.j)]_1
color.res[level.j & (res.i<res.j)]_3

color.w_rep(2,length(ydata))
color.w[level.i & (weight.i>weight.j)]_1
color.w[level.j & (weight.i<weight.j)]_3

if(isol==jsol)
{
plot(weight.i,col=color,xlab="Observations",ylab="Weights",main=paste("Weights of the root: ",isol))
}
else
{
if(isol>jsol)
{
plot(res.i,res.j,col=color.res,xlab=paste("Residuals of the root: ",isol),ylab=paste("Residuals of the root: ",jsol),main="Residuals")
abline(0,1)
}
else
{
plot(weight.i,weight.j,col=color.w,xlab=paste("Weights of the root: ",isol),ylab=paste("Weights of the root: ",jsol),main="Weights")
abline(0,1)
}
}
}
}


for(isol in 1:tot.sol)
{
x11()
par(mfcol=c(2,2))
y.fit.temp_y.fit[isol,]
res.temp_res[isol,]
weight.temp_weight[isol,]

level_weight.temp>=level.weight
color_rep(2,length(ydata))
color[level]_1

plot(y.fit.temp,res.temp,col=color,xlab="Fitted values",ylab="Residuals")

plot(y.fit.temp,res.temp*weight.temp,col=color,xlab="Fitted values",ylab="Weighted residuals")

qqnorm(res.temp,col=color)
qqline(res.temp)

qqnorm(res.temp*weight.temp,col=color)
qqline(res.temp*weight.temp)

##plot(weight.temp,col=color)
}

}
else
{

x11()
par(mfcol=c(1,1))

level.i_weight>=level.weight

color_rep(2,length(ydata))
color[level.i]_3

plot(weight,col=color,xlab="Observations",ylab="Weights",main=paste("Weights of the root"))

x11()
par(mfcol=c(2,2))

level_weight>=level.weight
color_rep(2,length(ydata))
color[level]_1

plot(y.fit,res,col=color,xlab="Fitted values",ylab="Residuals")

plot(y.fit,res*weight,col=color,xlab="Fitted values",ylab="Weighted residuals")

qqnorm(res,col=color)
qqline(res)

qqnorm(res*weight,col=color)
qqline(res*weight)

##plot(weight.temp,col=color)

}

}



