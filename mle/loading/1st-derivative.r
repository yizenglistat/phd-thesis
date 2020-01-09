"
This function computes first-order derivative of risk probability matrix with respect to beta (and delta)

Input:
  
  X       : features matrix
  beta    : feature coefficients
  delta   : gumbel parameter
  alpha1  : spline coefficients
  alpha2  : spline coefficients
  ord     : order of spline
  niknots : number of interior knots
  wrt     : derivative with respect to which variable, 'both' or only 'beta'

Output:
  first derivative matrix for beta (and delta), dimension is N x p or N x p+1
"
pmatrix_d<- function(X,beta,delta,alpha1,alpha2,ord,niknots,wrt='both'){
  # individual risk probability for disease type 1/type 2
  u1 <- g(X%*%beta,alpha1,ord,niknots)
  u2 <- g(X%*%beta,alpha2,ord,niknots)
  # first derivative of spline function (individual risk scores) for disease type 1/type 2    
  sp_d1 <- sp_d(X%*%beta,alpha1,ord,niknots)   
  sp_d2 <- sp_d(X%*%beta,alpha2,ord,niknots)
  # first derivative of probability of having two diseases simutaneously wrt beta
  p11_dbeta <- as.vector(-exp(-((-log(u1))^(1/delta)+(-log(u2))^(1/delta))^delta)*
                       ((-log(u1))^(1/delta)+(-log(u2))^(1/delta))^(delta-1)*
                       ((-log(u1))^(1/delta)/(log(u1))*(1-u1)*sp_d1
                         +(-log(u2))^(1/delta)/(log(u2))*(1-u2)*sp_d2))*X
  # first derivative of probability of having two diseases simutaneously wrt delta
  p11_ddelta<- -((-log(u2))^(1/delta)+(-log(u1))^(1/delta))^delta*
    ((-(-log(u2))^(1/delta)*log(-log(u2))/delta^2-(-log(u1))^(1/delta)*
        log(-log(u1))/delta^2)*delta/((-log(u2))^(1/delta)+(-log(u1))^(1/delta))+
       log((-log(u2))^(1/delta)+(-log(u1))^(1/delta)))*exp(-((-log(u2))^(1/delta)+(-log(u1))^(1/delta))^delta)
  # first derivative of risk probability matrix wrt beta and delta
  if(wrt=='both'){
    p11_d<- cbind(p11_dbeta,p11_ddelta)
    p10_d<- cbind(as.vector(u1*(1-u1)*sp_d1)*X, matrix(0,N,1))-p11_d
    p01_d<- cbind(as.vector(u2*(1-u2)*sp_d2)*X, matrix(0,N,1))-p11_d
    p00_d<- -p11_d-p10_d-p01_d
  }
  # first derivative of risk probability matrix wrt beta
  if(wrt=='beta'){
    p11_d<- p11_dbeta
    p10_d<- as.vector(u1*(1-u1)*sp_d1)*X-p11_d
    p01_d<- as.vector(u2*(1-u2)*sp_d2)*X-p11_d
    p00_d<- -p11_d-p10_d-p01_d
  }
  return(list(p00_d=p00_d,p10_d=p10_d,p01_d=p01_d,p11_d=p11_d))
}