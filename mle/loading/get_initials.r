# get initials
get_initials <- function(X,data,r,m){
  
  ini.beta <- rep(0,ncol(X)-1)#runif(ncol(X)-1,-1,1)
  #ini.beta <- ini.beta/sqrt(sum(ini.beta^2))
  #ini.beta[1]<-ifelse(ini.beta[1]>0,ini.beta[1],-ini.beta[1])
  
  u<-sort(X%*%c(1,ini.beta))
  #u<-sort(X%*%ini.beta)
  
  #alpha1 <- bayesglm(data[,1]~Bsp(u,r,m)-1,family = binomial())$coef
  #alpha2 <- bayesglm(data[,2]~Bsp(u,r,m)-1,family = binomial())$coef
  #alpha1= c(0.2471403,  0.0000001,  0.5081673 , 0.1195941 , 0.2373479,-4.5794987)
  #alpha2=c(0.50762103 , 0.03779975,  0.64315908,  0.09101801  ,0.54073814,-4.95008617)
  #ini.alpha <- c(alpha1,alpha2)
  #ini.alpha = runif(r+m+r+m)
  
  library(lsei)
  x = seq(-10,10,length.out=100)+rnorm(100,0,0.4)
  ini.alpha = pnnls(Bsp(x,r,m),x,1)$x
  ini.alpha=c(ini.alpha,ini.alpha)
  
  ini.delta <- runif(1,0.001,1-0.001)
  return(list(ini.beta=ini.beta,
              ini.alpha=ini.alpha,
              ini.delta=ini.delta))
}


