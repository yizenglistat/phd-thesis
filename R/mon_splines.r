Bsp <- function(u,r,m){
  xi=quantile(u,probs = seq(0,1,length.out = m-1))
  xi[1] = xi[1]+0.001
  xi[length(xi)]=xi[length(xi)]-0.001
  out<-iSpline(knots=xi, degree=r-1, x=u, intercept = T)
  out = cbind(1,out)
  return(out)
}

eta.Bsp <- function(u,alpha,r,m){
  return(Bsp(u,r,m)%*%alpha)
}

# inverse of logit function
logistic <- function(u){return(1/(1+exp(-u)))}

# generilized link function, g_k in our case, depend on alpha_k
g <- function(u,alphak,r,m){return(logistic(eta.Bsp(u,alphak,r,m)))}

# Gumbel copula
gumbel <- function(u,v,delta){
  return(exp(-((-log(u))^(1/delta)+(-log(v))^(1/delta))^(delta)))
}

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

