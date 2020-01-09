#***********************************************************************#
#                        GEM algorithm function  
#***********************************************************************#
source('./loading/constrOptimNL.r')
MLE<- function(X, cj, data, Se, Sp, r, m, 
               eps, verbose=T, maxiter=20){
  # setting
  iter <- 1
  ini.err <- 1
  tem.err <- eps+1
  ctime <- proc.time()
  cost <- c()
  # initials
  initials <- get_initials(X,data,r,m)
  ini.alpha <- initials$ini.alpha
  ini.beta <- initials$ini.beta
  ini.delta <- initials$ini.delta
  
  # size of parameters
  nalpha <- length(ini.alpha)
  nbeta <- length(ini.beta)
  ndelta <- length(ini.delta)
  # EM algorithm iteration
  while(iter < maxiter){

    # inital parameters
    ini.par <- c(ini.alpha,ini.beta,ini.delta) # c(ini.alpha,ini.beta,ini.delta)
    # E-step 
    pp <- p.matrix(X,ini.alpha,ini.beta,ini.delta,r,m)
    tau.Y <- tau.Y1Y2(data,pp,Se,Sp,cj)$tau
    ini.obj <- E.loglik(ini.alpha,ini.delta,ini.beta,X,data,r,m,tau.Y)
    cost <- c(cost,ini.obj)
    # M-step
    
    # transform to non-negaitve argument
    #ini.beta[1] <- log(ini.beta[1])
    # rep=1
    # nest_err=1
    # while(nest_err>1e-1){
    pp <- p.matrix(X,ini.alpha,ini.beta,ini.delta,r,m)
    tau.Y <- tau.Y1Y2(data,pp,Se,Sp,cj)$tau
    
    tem.beta <- optim(par=ini.beta, fn=l.beta,
                      tem.alpha=ini.alpha,tem.delta=ini.delta,
                      X=X,data=data,r=r,m=m,tau.Y=tau.Y)$par
    # tem.beta[1] <- exp(tem.beta[1])
    w = 1.5
    tem.beta = (1+w)*tem.beta - w*ini.beta
    
    #tem.beta = ini.beta + w * (tem.beta-ini.beta)
    #tem.beta <- tem.beta / sqrt(sum(tem.beta^2))
    print("true beta")
    print(beta)
    print("beta updated")
    print(tem.beta)
    
    tem.delta <- nlminb(start=ini.delta, objective=l.delta, 
                       lower = 0.001,upper=0.999,
                       tem.alpha=ini.alpha,ini.beta=tem.beta,tau.Y=tau.Y, 
                       X=X,data=data,r=r,m=m)$par
    print("delta updated")
    print(tem.delta)
    #nest_err = sqrt(sum((c(ini.beta,ini.delta)-c(tem.beta,tem.delta))^2))
    #print(nest_err)
    ini.beta = tem.beta
    ini.delta = tem.delta
  #}
  rep = 1
  for(rep in 1:2){
    pp <- p.matrix(X,ini.alpha,ini.beta,ini.delta,r,m)
    tau.Y <- tau.Y1Y2(data,pp,Se,Sp,cj)$tau
    # log-transform
    ini.alpha[2:(r+m)]<-log(ini.alpha[2:(r+m)])
    ini.alpha[(r+m+2):(r+m+r+m)]<-log(ini.alpha[(r+m+2):(r+m+r+m)])
    
    tem.alpha <- optim(par=ini.alpha, fn=l.alpha,
                        ini.beta=ini.beta,ini.delta=ini.delta,
                        X=X,data=data,r=r,m=m,tau.Y=tau.Y)$par
    
    #tem.alpha = ini.alpha + w*(tem.alpha-ini.alpha)
    w = 1.5
    tem.alpha = (1+w)*tem.alpha-w*ini.alpha
    
    # transform to non-negative parameters
    ini.alpha[2:(r+m)]<-exp(ini.alpha[2:(r+m)])
    ini.alpha[(r+m+2):(r+m+r+m)]<-exp(ini.alpha[(r+m+2):(r+m+r+m)])
    
    tem.alpha[2:(r+m)]<-exp(tem.alpha[2:(r+m)])
    tem.alpha[(r+m+2):(r+m+r+m)]<-exp(tem.alpha[(r+m+2):(r+m+r+m)])
  
    # show the iteration details and plots
    if (verbose){
      print(paste0(paste0(rep('*',30),collapse =''),iter,"th Iteration",paste0(rep('*',30),collapse ='')))
      print("beta")
      print(tem.beta)
      print("alpha1")
      print(as.vector(tem.alpha)[1:(nalpha/2)])
      print("alpha2")
      print(as.vector(tem.alpha)[(nalpha/2+1):nalpha])
      print("delta")
      print(tem.delta)
      print("err")
      print(tem.err)
      print("cost")
      print(cost)
      # generate the interation plot
      if (iter==1){
        par(mfrow=c(1,2))
        u<-sort(X%*%c(1,ini.beta))
        #u<-sort(X%*%tem.beta)
        plot(u,g(u,tem.alpha[1:(nalpha/2)],r,m),type='l',col='gray',
             ylab=expression(hat(g)[1](u)))
        plot(u,g(u,tem.alpha[(nalpha/2+1):nalpha],r,m),type='l',col='gray',
             ylab=expression(hat(g)[2](u)))
      }else{
        par(mfrow=c(1,2))
        u<-sort(X%*%c(1,ini.beta))
        #u<-sort(X%*%ini.beta)
        plot(u,g(u,ini.alpha[1:(nalpha/2)],r,m),type='l',col='gray',ylab=expression(hat(g)[1](u)))
        lines(u,g(u,tem.alpha[1:(nalpha/2)],r,m),type='l',col='black')
        lines(u,links[[1]](u),type='l',col='green')
        plot(u,g(u,ini.alpha[(nalpha/2+1):nalpha],r,m),type='l',col='gray',ylab=expression(hat(g)[2](u)))
        lines(u,g(u,tem.alpha[(nalpha/2+1):nalpha],r,m),type='l',col='black')
        lines(u,links[[2]](u),type='l',col='green')
      }
    }
    
    print("alpha updated")
    ini.alpha = tem.alpha
  }
    
    

    
    # update the error
    tem.par <- c(tem.alpha,tem.beta,tem.delta) # c(tem.alpha,tem.beta,tem.delta)
    pp <- p.matrix(X,tem.alpha,tem.beta,tem.delta,r,m)
    tau.Y <- tau.Y1Y2(data,pp,Se,Sp,cj)$tau
    tem.obj <- E.loglik(tem.alpha,tem.delta,tem.beta,X,data,r,m,tau.Y)
    tem.err <- abs(tem.obj-ini.obj)  # sqrt(sum((tem.par-ini.par)^2))
    
    # check the convergence
    #if (tem.err<eps){break}
    
    #if (abs(ini.err-tem.err)<1e-3){break}
    
    # update the loop
    iter <- iter + 1 
    ini.alpha <- tem.alpha
    ini.beta <- tem.beta
    ini.delta <- tem.delta
    ini.err <- tem.err
  }
  
  # estimate the execution time
  time <- paste0(round(unname(((proc.time() - ctime)/60)[1]),3)," mins")
  print(paste0(paste0(rep('*',30),collapse =''),"Output Below",paste0(rep('*',30),collapse ='')))
  # return the estimates for paramters
  return(list(iterations=iter,
              cost=cost,
              time=time,
              final.alpha1=as.vector(tem.alpha)[1:(nalpha/2)],
              final.alpha2=as.vector(tem.alpha)[(nalpha/2+1):nalpha],
              final.beta=tem.par[1:nbeta],
              final.delta=tail(as.vector(tem.par),n=1)
  ))
}