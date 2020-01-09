source('./loading/optim_alpha_beta.r')
#***********************************************************************#
#                   Optimization alpha function
#***********************************************************************#
# function for beta optimization
l.alpha<- function(tem.par,ini.delta,ini.beta,
                   X,data,r,m,tau.Y){
  # transform to nonnegative parameters
  
  tem.par[2:(r+m)]<-exp(tem.par[2:(r+m)])
  tem.par[(r+m+2):(r+m+r+m)]<-exp(tem.par[(r+m+2):(r+m+r+m)])
  
  
  out <- -E.loglik(tem.par,ini.delta,ini.beta,
                   X,data,r,m,tau.Y)
  return(out)
  #penalty <- 1
  # R <- matrix(NA,r+m-2,r+m)
  # for (row_ind in 1:(r+m-2)){
  #   R[row_ind, ] <- c(rep(0,row_ind-1),1,-2,1,rep(0,r+m-row_ind-2)) 
  # }
  # 
  # penalty_fun <- function(u){
  #   return(t(Bsp_hess(u,r,m))%*%Bsp_hess(u,r,m))
  # }
  
  
  
  # #M <- penalty_fun(X%*%ini.beta)
  # a<-range(X%*%ini.beta)[1]
  # b<-range(X%*%ini.beta)[2]
  # n<-100
  # h<-(b-a)/n
  # u<-seq(a, b, by=h)
  # y<-Bsp_hess(u,r,m)
  # omg<-list()
  # mat<-matrix(NA,r+m,r+m)
  # 
  # for (i in 1:dim(y)[1]){
  #   for (j in 1:dim(y)[2]){
  #     for (k in 1:dim(y)[2]){
  #      mat[j,k]<-y[i,j]*y[i,k]
  #     }
  #   }
  #   omg[[i]]<-mat
  # }
  # 
  # ss<-matrix(0,r+m,r+m)
  # for(i in 2:n){
  #   ss<-ss+omg[[i]]
  # }
  
  #pen_mat <- h*(omg[[1]]/2+ss+omg[[n+1]]/2)
  #pen1 <- t(tem.par[1:(r+m)])%*%pen_mat%*%tem.par[1:(r+m)] 
  #pen2 <- t(tem.par[(r+m+1):(r+m+r+m)])%*%pen_mat%*%tem.par[(r+m+1):(r+m+r+m)]
  #out <- out + 1 * sum((R%*%tem.par[1:(r+m)])^2) + 2 * sum((R%*%tem.par[(r+m+1):(r+m+r+m)])^2)
  
  
  #out <- out + 10 * penalty * sum((R%*%tem.par[1:(r+m)])^2) 
  #           + 100 * penalty * sum((R%*%tem.par[(r+m+1):(r+m+r+m)])^2)
  #out <- out #+ penalty * sqrt(sum((Bsp_hess(X%*%ini.beta,r,m)%*%tem.par[1:(r+m)])^2))
            #+ 10 * penalty * sqrt(sum((Bsp_hess(X%*%ini.beta,r,m)%*%tem.par[(r+m+1):(r+m+r+m)])^2)) 
  return(out)
}
