"
This function computes second-order derivative of risk probability matrix with respect to beta (and delta)

Input:
  
  X       : features matrix
  beta    : feature coefficients
  delta   : gumbel parameter
  alpha1  : spline coefficients
  alpha2  : spline coefficients
  ord     : order of spline
  niknots : number of interior knots
  wrt     : derivative with respect to which variable, 'both' or just 'beta'

Output:
  second derivative matrix for beta (and delta), dimension is N x (p x p) or N x (p+1 x p+1)
"
pp.dd<- function(X,beta,delta,alpha1,alpha2,ord,niknots,wrt='both'){
  N <- nrow(X)
  # individual risk probability for disease type 1/type 2
  u1 <- g(X%*%beta,alpha1,ord,niknots)
  u2 <- g(X%*%beta,alpha2,ord,niknots)
  # first derivative of spline function (individual risk scores) for disease type 1/type 2    
  sp_d1 <- sp_d(X%*%beta,alpha1,ord,niknots)   
  sp_d2 <- sp_d(X%*%beta,alpha2,ord,niknots)
  # second derivative of spline function (individual risk scores) for disease type 1/type 2    
  sp_dd1 <- sp_dd(X%*%beta,alpha1,ord,niknots)   
  sp_dd2 <- sp_dd(X%*%beta,alpha2,ord,niknots)
  # split X matrix into a list and each element is a row (subject) in the list
  X_list  <- lapply(seq_len(nrow(X)), function(idx) X[idx,])
  # a list of length = N, each element is a matrix of (p x p), i.e., X[i,]%*%t(X[i,])
  Xsq_list <- lapply(X_list,function(x){(x%*%t(x))})
  # second derivative w.r.t beta 
  p11_ddbeta <- exp(-((-log(u1))^(1/delta)+(-log(u2))^(1/delta))^delta)*
        ((-log(u1))^(1/delta)+(-log(u2))^(1/delta))^(delta-2)*
        (
          ((-log(u1))^(1/delta)+(-log(u2))^(1/delta))^(delta)*
          ((-log(u1))^(1/delta-1)*(1-u1)*sp_d1+(-log(u2))^(1/delta-1)*(1-u2)*sp_d2)^2+
          ((delta-1)/delta)*((-log(u1))^(1/delta-1)*(1-u1)*sp_d1+(-log(u2))^(1/delta-1)*(1-u2)*sp_d2)^2+
          ((-log(u1))^(1/delta)+(-log(u2))^(1/delta))*
          (
            (1/delta-1)*(-log(u1))^(1/delta-2)*(1-u1)^2*sp_d1^2+
            (-log(u1))^(1/delta-1)*u1*(1-u1)*sp_d1^2-
            (-log(u1))^(1/delta-1)*(1-u1)*sp_dd1+
            (1/delta-1)*(-log(u2))^(1/delta-2)*(1-u2)^2*sp_d2^2+
            (-log(u2))^(1/delta-1)*u2*(1-u2)*sp_d2^2-
            (-log(u2))^(1/delta-1)*(1-u2)*sp_dd2
          )
        )
  p11_ddbeta <- listOps(p11_ddbeta,Xsq_list,"*")
  # second derivative w.r.t beta and delta
  p11_ddbetadelta <- 
    as.vector(exp(-((-log(u2))^(1/delta)+(-log(u1))^(1/delta))^delta)*
    ((-log(u2))^(1/delta)+(-log(u1))^(1/delta))^delta*
    ((-(-log(u2))^(1/delta)*log(-log(u2))/delta^2-(-log(u1))^(1/delta)*log(-log(u1))/delta^2)*delta/((-log(u2))^(1/delta)+(-log(u1))^(1/delta))+
      log((-log(u2))^(1/delta)+(-log(u1))^(1/delta)))*
    ((-log(u1))^(1/delta)+(-log(u2))^(1/delta))^(delta-1)*
    ((-log(u1))^(1/delta)/(log(u1))*(1-u1)*sp_d1+(-log(u2))^(1/delta)/(log(u2))*(1-u2)*sp_d2))*X+
    as.vector(-exp(-((-log(u1))^(1/delta)+(-log(u2))^(1/delta))^delta)*
    ((-log(u2))^(1/delta)+(-log(u1))^(1/delta))^(delta-1)*
    ((-(-log(u2))^(1/delta)*log(-log(u2))/delta^2-(-log(u1))^(1/delta)*log(-log(u1))/delta^2)*(delta-1)/((-log(u2))^(1/delta)+(-log(u1))^(1/delta))+
      log((-log(u2))^(1/delta)+(-log(u1))^(1/delta)))*
    ((-log(u1))^(1/delta)/(log(u1))*(1-u1)*sp_d1+(-log(u2))^(1/delta)/(log(u2))*(1-u2)*sp_d2))*X+
    as.vector(-exp(-((-log(u1))^(1/delta)+(-log(u2))^(1/delta))^delta)*
    ((-log(u1))^(1/delta)+(-log(u2))^(1/delta))^(delta-1)*
    ((-log(u1))^(1/delta-1)*log(-log(u1))*(1/delta^2)*(1-u1)*sp_d1+(-log(u2))^(1/delta-1)*log(-log(u2))*(1/delta^2)*(1-u2)*sp_d2))*X
  p11_ddbetadelta <- lapply(1:N,function(x){matrix(p11_ddbetadelta[x,],ncol(X),1)})
  # second derivative w.r.t delta
  p11_dddelta <- 
    exp(-((-log(u2))^(1/delta)+(-log(u1))^(1/delta))^delta)*
    ((-log(u2))^(1/delta)+(-log(u1))^(1/delta))^(2*delta)*
    ((-(-log(u2))^(1/delta)*log(-log(u2))/delta^2-(-log(u1))^(1/delta)*
          log(-log(u1))/delta^2)*delta/((-log(u2))^(1/delta)+(-log(u1))^(1/delta))+
         log((-log(u2))^(1/delta)+(-log(u1))^(1/delta)))^2-
    exp(-((-log(u2))^(1/delta)+(-log(u1))^(1/delta))^delta)*
    ((-log(u2))^(1/delta)+(-log(u1))^(1/delta))^delta*
    ((-(-log(u2))^(1/delta)*log(-log(u2))/delta^2-(-log(u1))^(1/delta)*
          log(-log(u1))/delta^2)*delta/((-log(u2))^(1/delta)+(-log(u1))^(1/delta))+
         log((-log(u2))^(1/delta)+(-log(u1))^(1/delta)))^2-
    exp(-((-log(u2))^(1/delta)+(-log(u1))^(1/delta))^delta)*
    ((-log(u2))^(1/delta)+(-log(u1))^(1/delta))^delta*
    (
      2*(-(-log(u1))^(1/delta)*log(-log(u1))/delta^2-(-log(u2))^(1/delta)*log(-log(u2))/delta^2)/((-log(u1))^(1/delta)+(-log(u2))^(1/delta))-
      delta*(-(-log(u1))^(1/delta)*log(-log(u1))/delta^2-(-log(u2))^(1/delta)*log(-log(u2))/delta^2)^2/((-log(u1))^(1/delta)+(-log(u2))^(1/delta))^2+
      delta*(-(-log(u1))^(1/delta)*(log(-log(u1)))^2/delta^4-(-log(u2))^(1/delta)*(log(-log(u2)))^2/delta^4)/((-log(u1))^(1/delta)+(-log(u2))^(1/delta))+
      delta*(-2*(-log(u1))^(1/delta)*log(-log(u1))/delta^3-2*(-log(u2))^(1/delta)*log(-log(u2))/delta^3)/((-log(u1))^(1/delta)+(-log(u2))^(1/delta))
    )
  p11_dddelta <- lapply(1:N,function(x){matrix(p11_dddelta[x,],1,1)})
  # second derivative of risk probability matrix wrt beta and delta
  if(wrt=='both'){
    p11.dd <- mapply(function(p11_ddbeta,p11_ddbetadelta,p11_dddelta){rbind(cbind(p11_ddbeta,p11_ddbetadelta),cbind(t(p11_ddbetadelta),p11_dddelta))},
                     p11_ddbeta,p11_ddbetadelta,p11_dddelta,SIMPLIFY = FALSE)

    p10.dd <- listOps(u1*(1-u1)*((1-2*u1)*sp_d1^2+sp_dd1),Xsq_list,"*")
    p10.dd <- lapply(p10.dd,function(x){rbind(cbind(x,matrix(0,p,1)),matrix(0,1,p+1))})
    p10.dd <- listOps(p10.dd,p11.dd,"-")
    
    p01.dd <- listOps(u1*(1-u2)*((1-2*u2)*sp_d2^2+sp_dd2),Xsq_list,"*")
    p01.dd <- lapply(p01.dd,function(x){rbind(cbind(x,matrix(0,p,1)),matrix(0,1,p+1))})
    p01.dd <- listOps(p01.dd,p11.dd,"-")

    p00.dd <- mapply(function(x,y,z){-x-y-z},p10.dd,p01.dd,p11.dd,SIMPLIFY = FALSE)
  }
  # second derivative of risk probability matrix wrt beta
  if(wrt=='beta'){
    p11.dd <- p11_ddbeta
    p10.dd <- listOps(u1*(1-u1)*((1-2*u1)*sp_d1^2+sp_dd1),Xsq_list,"*")
    p10.dd <- listOps(p10.dd,p11.dd,"-")
    
    p01.dd <- listOps(u1*(1-u2)*((1-2*u2)*sp_d2^2+sp_dd2),Xsq_list,"*")
    p01.dd <- listOps(p01.dd,p11.dd,"-")

    p00.dd <- mapply(function(x,y,z){-x-y-z},p10.dd,p01.dd,p11.dd,SIMPLIFY = FALSE)
  }
  return(list(p00.dd=p00.dd,p10.dd=p10.dd,p01.dd=p01.dd,p11.dd=p11.dd))
}