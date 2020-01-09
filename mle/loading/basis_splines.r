#***********************************************************************#
#                       B-splines functions
#***********************************************************************#
# eta.Bsp(u,alpha,r,m) single-index risk score for kth infection, depend on alpha
# alpha determine the kind of infection
# r: the order of B-splines, i.e., degree is p-1, e.g. r = 3 is a quadratic Bspline
# m: numerber of interior knots
library(survival)
#i-spline basis functions evaluating at u
Bsp <- function(u,r,m){
  # normalize the data in order to fix knots
  #x<-scale(u)
  #F0<-pnorm(x)
  #xi<-(1:(m-1))/m;
  #xi <- c(0.05,0.1,0.15,0.2,0.25,0.5,0.7,0.8,0.9)
  #xi<-c( rep(min(xi),r-1), xi, rep(max(xi),r-1) )
  xi=quantile(u,probs = seq(0,1,length.out = m-1))
  xi[1] = xi[1]+0.001
  xi[length(xi)]=xi[length(xi)]-0.001
  #print(xi)
  #print(xi)
  out<-iSpline(knots=xi, degree=r-1, x=u, intercept = T)
  # add a constant function
  out = cbind(1,out)
  return(out)
}

eta.Bsp <- function(u,alpha,r,m){
  return(Bsp(u,r,m)%*%alpha)
}


