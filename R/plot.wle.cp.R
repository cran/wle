#############################################################
#                                                           #
#	PLOT.WLE.CP function                                #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: October, 10, 1999                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 1999 Claudio Agostinelli              #
#                                                           #
#############################################################

plot.wle.cp_function(object,base.line=0,num.max=20,plot.it=TRUE,log.scale=FALSE)
{
wcp_object$wcp
num.model_dim(wcp)[1]
if(num.model<num.max)
{
warning("plot.wle.cp: The number of model is less than num.max")
num.max_num.model
}
if(num.max<2){num.max_2}

nvar_dim(wcp)[2]-1
good.model_(apply(wcp[,1:nvar],1,sum)+base.line>=wcp[,nvar+1])
wcp.good_matrix(wcp[good.model,],ncol=nvar+1)
wcp.bad_matrix(wcp[!good.model,],ncol=nvar+1)
ordine.good_order(wcp.good[,nvar+1])
ordine.bad_order(wcp.bad[,nvar+1])
wcp.good_matrix(wcp.good[ordine.good,],ncol=(nvar+1))
num.good_dim(wcp.good)[1]
wcp.bad_matrix(wcp.bad[ordine.bad,],ncol=(nvar+1))
num.bad_dim(wcp.bad)[1]

label.good_character()
for(i in 1:nvar){
label.good_paste(label.good,wcp.good[,i],sep="")
}
label.bad_character()
for(i in 1:nvar){
label.bad_paste(label.bad,wcp.bad[,i],sep="")
}


xcoord.good_apply(matrix(wcp.good[,1:nvar],ncol=nvar),1,sum)[1:min(num.max,num.good)]
ycoord.good_wcp.good[,nvar+1][1:min(num.max,num.good)]

label.good_label.good[1:min(num.max,num.good)]

xcoord.best_xcoord.good[1]
ycoord.best_ycoord.good[1]

label.best_label.good[1]

if(length(xcoord.good)==1)
{
xcoord.good_0
ycoord.good_0
plot.good_FALSE
}else
{
xcoord.good_xcoord.good[-1]
ycoord.good_ycoord.good[-1]
label.good_label.good[-1]
plot.good_TRUE
}

if(num.max>num.good)
{
xcoord.bad_apply(matrix(wcp.bad[,1:nvar],ncol=nvar),1,sum)[1:min(num.bad,num.max-num.good)]
ycoord.bad_wcp.bad[,nvar+1][1:min(num.bad,num.max-num.good)]
label.bad_label.bad[1:min(num.bad,num.max-num.good)]
plot.bad_TRUE
}else
{
xcoord.bad_0
ycoord.bad_0
plot.bad_FALSE
}

xlim.min_min(xcoord.good,xcoord.bad,xcoord.best)
xlim.max_max(xcoord.good,xcoord.bad,xcoord.best)

yetichetta_"WCp"

if(log.scale)
{
ycoord.good_log10(ycoord.good+min(ycoord.good,ycoord.bad,ycoord.best)+1)
ycoord.bad_log10(ycoord.bad+min(ycoord.good,ycoord.bad,ycoord.best)+1)
ycoord.best_log10(ycoord.best+min(ycoord.good,ycoord.bad,ycoord.best)+1)
yetichetta_"WCp log10 scale"
}

ylim.min_min(ycoord.good,ycoord.bad,ycoord.best)
ylim.max_max(ycoord.good,ycoord.bad,ycoord.best)


if(plot.it)
{
plot(xcoord.best,ycoord.best,xlim=c(xlim.min,xlim.max),ylim=c(ylim.min,ylim.max),xlab="Number of Predictors",ylab=yetichetta,type="n")
text(xcoord.best,ycoord.best,col=4,labels=label.best)

if(plot.good)
{
text(xcoord.good,ycoord.good,col=3,labels=label.good)
}

if(plot.bad)
{
text(xcoord.bad,ycoord.bad,col=2,labels=label.bad)
}


if(!log.scale)
{
abline(base.line,1,col=2)
abline(0,1)
}
else
{
vettx_seq(xlim.min,xlim.max,0.5)
vetty_log10(vettx+min(ycoord.good,ycoord.bad,ycoord.best)+1)
vetty.base.line_log10(vettx+min(ycoord.good,ycoord.bad,ycoord.best)+1+base.line)
lines(vettx,vetty.base.line,col=2,type="l")
lines(vettx,vetty,type="l")
}


}

invisible(list(num.good=num.good,num.bad=num.bad,wcp.good=wcp.good, wcp.bad=wcp.bad))
}
