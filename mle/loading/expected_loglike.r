# expectation for log-likelihood function, the objective function
E.loglik<- function(obj.alpha,obj.delta,obj.beta,X,data,r,m,tau.Y){
  # Q1 expectation of log-likelihood function
  penalty <- -1e20
  
  tem.p <- p.matrix(X,obj.alpha,obj.beta,obj.delta,r,m)

  if (any(tem.p<0) | any(is.na(tem.p)) | any(is.nan(tem.p))){
    return(penalty)
  }
  
  
  
  # objetive function
  E.loglik.Q1 <- sum(tau.Y*(log(tem.p)))
  
  return(E.loglik.Q1)
}