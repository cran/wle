#############################################################
#                                                           #
#	PLOT.WLE.LM function                                #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: December, 19, 2000                            #
#	Version: 0.3                                        #
#                                                           #
#	Copyright (C) 2000 Claudio Agostinelli              #
#                                                           #
#############################################################

plot.wle.lm <- function(object, level.weight=0.5, ask = interactive() && .Device != "postscript") {

old.par <- par(no.readonly=TRUE)
on.exit(par(old.par))

if (!inherits(object, "wle.lm")) stop("Use only with 'wle.lm' objects")

if (ask) par(ask = TRUE)
 
if (!is.null(object$model)) {
ydata <- object$model[,1]
xdata <- object$model[,-1]
} else if (!is.null(object$x) & !is.null(object$y)) {
ydata <- object$y
xdata <- object$x
} else {
stop("Please, rerun wle.lm with model=TRUE (or x=TRUE and y=TRUE)")
}

if (level.weight<0 | level.weight >1) {
cat("plot.wle.lm: level.weight should be between zero and one, set to 0.5 \n")
level.weight <- 0.5
}

param <- as.matrix(object$coefficients)
res <- as.matrix(object$residuals)
y.fit <- as.matrix(object$fitted.values)
weight <- as.matrix(object$weights)
tot.weight <- object$tot.weights
tot.sol <- object$tot.sol

if (tot.sol>1) {
par(mfcol=c(tot.sol,tot.sol))

for (isol in 1:tot.sol) {
for (jsol in 1:tot.sol) {
y.fit.i <- y.fit[isol,]
res.i <- res[isol,]
weight.i <- weight[isol,]

y.fit.j <- y.fit[jsol,]
res.j <- res[jsol,]
weight.j <- weight[jsol,]

level.i <- weight.i>=level.weight
level.j <- weight.j>=level.weight

color <- rep(2,length(ydata))
color[level.i] <- 1
color[level.j] <- 3

color.res <- rep(2,length(ydata))
color.res[level.i & (res.i>res.j)] <- 1
color.res[level.j & (res.i<res.j)] <- 3

color.w <- rep(2,length(ydata))
color.w[level.i & (weight.i>weight.j)] <- 1
color.w[level.j & (weight.i<weight.j)] <- 3

if (isol==jsol) {
plot(weight.i,col=color,xlab="Observations",ylab="Weights",main=paste("Weights of the root: ",isol))
} else {
if (isol>jsol) {
plot(res.i,res.j,col=color.res,xlab=paste("Residuals of the root: ",isol),ylab=paste("Residuals of the root: ",jsol),main="Residuals")
abline(0,1)
} else {
plot(weight.i,weight.j,col=color.w,xlab=paste("Weights of the root: ",isol),ylab=paste("Weights of the root: ",jsol),main="Weights")
abline(0,1)
}
}
}
}


for (isol in 1:tot.sol) {
par(mfcol=c(2,2))
y.fit.temp <- y.fit[isol,]
res.temp <- res[isol,]
weight.temp <- weight[isol,]

level <- weight.temp>=level.weight
color <- rep(2,length(ydata))
color[level] <- 1

plot(y.fit.temp,res.temp,col=color,xlab="Fitted values",ylab="Residuals")

plot(y.fit.temp,res.temp*weight.temp,col=color,xlab="Fitted values",ylab="Weighted residuals")

qqnorm(res.temp,col=color)
qqline(res.temp)

qqnorm(res.temp*weight.temp,col=color)
qqline(res.temp*weight.temp)
}
} else {

par(mfcol=c(1,1))
level.i <- weight>=level.weight
color <- rep(2,length(ydata))
color[level.i] <- 3

plot(weight,col=color,xlab="Observations",ylab="Weights",main="Weights of the root")

par(mfcol=c(2,2))

level <- weight>=level.weight
color <- rep(2,length(ydata))
color[level] <- 1

plot(y.fit,res,col=color,xlab="Fitted values",ylab="Residuals")

plot(y.fit,res*weight,col=color,xlab="Fitted values",ylab="Weighted residuals")

qqnorm(res,col=color)
qqline(res)

qqnorm(res*weight,col=color)
qqline(res*weight)
}
}



