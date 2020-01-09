#***********************************************************************#
#                   Optimization delta function
#***********************************************************************#
# function for delta optimization
l.delta<- function(tem.par,tem.alpha,ini.beta,
                   X,data,r,m,tau.Y){
  return(-E.loglik(tem.alpha,tem.par,ini.beta,
                   X,data,r,m,tau.Y))
}