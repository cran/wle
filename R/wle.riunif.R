#############################################################
#                                                           #
#	wle.riunif function                                 #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: January, 19, 2001                             #
#	Version: 0.1                                        #
#                                                           #
#############################################################

wle.riunif <- function(x,inf,sup) {
    round(runif(x,inf-0.4999,sup+0.4999))
}






