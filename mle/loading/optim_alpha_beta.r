l.alpha_beta<- function(tem.par,ini.delta,
                   X,data,r,m,tau.Y){
  # transform to nonnegative parameters
  tem.par[1:(r+m-1)]<-exp(tem.par[1:(r+m-1)])
  tem.par[(r+m+1):(r+m+r+m-1)]<-exp(tem.par[(r+m+1):(r+m+r+m-1)])
  
  myalpha = tem.par[1:(r+m+r+m)]
  mybeta = tem.par[(r+m+r+m+1):(r+m+r+m+ncol(X))]
  
  out <- -E.loglik(myalpha,ini.delta,mybeta,
                   X,data,r,m,tau.Y)
  return(out)
}