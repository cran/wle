#############################################################
#                                                           #
#	WLE.SMOOTH function                                 #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: October, 10, 1999                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 1999 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.smooth_function(weight=0.31,costant=3,level=0.2,dimension=1,raf=1,interval=c(0.00001,0.5),tol=10^-6,max.iter=1000)
{

delta_function(smooth,costant,level,dimension){
level*(((smooth+1)/smooth)^(dimension/2)*exp(costant^2/(2*(dimension+1)))-1)}

if(raf==3) {w_function(smooth,costant,level,dimension,weight){(1-((delta(smooth,costant,level,dimension)**2)/((delta(smooth,costant,level,dimension)**2) + 2)))-weight}
}
else
{
if(raf==2) adelta_function(d){2-(2+d)*exp(-d)}
else
{
if(raf==1) adelta_function(d){2*(sqrt(d+1)-1)}

}
w_function(smooth,costant,level,dimension,weight){
(adelta(delta(smooth,costant,level,dimension))+1)/(delta(smooth,costant,level,dimension)+1)-weight
}
}

return(uniroot(w,interval=interval,costant=costant,level=level,dimension=dimension,weight=weight,maxiter=max.iter,tol=tol))
}
