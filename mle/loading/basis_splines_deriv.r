# # first derivative of B-spline function of order r, evaluating at u
# eta.Bsp.deriv <- function(u,alpha,r,m){
#   # normalize the data in order to fix knots
#   x<-scale(u)
#   F0<-pnorm(x);
#   xi<-min(F0) + ((max(F0) - min(F0) + 0.0001)*(0:(m+1))/(m+1));
#   xi<-c( rep(min(xi),r-1), xi, rep(max(xi),r-1) )
#   alpha.deriv<- (r-1)*(diff(alpha))/(diff(xi,lag=r-1)[-c(1,m+r+1)])
#   return(Bsp(u,r-1,m)%*%alpha.deriv)
# }


# first derivative

sp_d <- function(u,alpha,ord,niknots){
  # normalize the data in order to fix knots
  x<-scale(u)
  F0<-pnorm(x);
  #F0 = u
 # xi<-min(F0) + ((max(F0) - min(F0) + 0.001)*(1:(niknots-1))/(niknots));
  #xi<-c( rep(min(xi),r-1), xi, rep(max(xi),r-1) )
  #out<-iSpline(knots=xi, degree=ord-1, x=F0, intercept=T,derivs = 1)
  xi = (1:(niknots-1))/(niknots);
  out = mSpline(knots=xi, degree=ord-1, x=F0, intercept=T)
  return(out%*%alpha[-length(alpha)])
}


# second derivative

sp_dd <- function(u,alpha,ord,niknots){
  # normalize the data in order to fix knots
  x<-scale(u)
  F0<-pnorm(x);
  #F0 = u
  #xi<-min(F0) + ((max(F0) - min(F0) + 0.001)*(1:(niknots-1))/(niknots));
  xi = (1:(niknots-1))/(niknots);
  #xi<-c( rep(min(xi),r-1), xi, rep(max(xi),r-1) )
  #out<-iSpline(knots=xi, degree=ord-1, x=F0, intercept=T,derivs = 1)
  out = mSpline(knots=xi, degree=ord-1, x=F0, intercept=T,derivs = 1)
  return(out%*%alpha[-length(alpha)])
}