#***********************************************************************#
#                    Optimization beta function
#***********************************************************************#
# function for beta optimization
l.beta<- function(tem.par,tem.alpha,tem.delta,
                  X,data,r,m,tau.Y){
  #tem.par[1]<-exp(tem.par[1])
  return(-E.loglik(tem.alpha,tem.delta,tem.par,
                   X,data,r,m,tau.Y))
}