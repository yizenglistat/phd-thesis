"
Maximum Likelihood Estimation via adjusted EM algorithm

Args:

  X       : feature matrix
  cj      : designed group size vector
  data    : observed response matrix
  Se      : sensitivity vector
  Sp      : specificity vector
  r       : order of splines
  m       : number of interior knots
  verbose : TRUE or FALSE, presents iteration plots and details during estimation

Return:

  A list of estimates and costs iteration history

"
MLE <- function(X, cj, data, Se, Sp, r, m, verbose=TRUE){
  
  # initial setting  
  initials <- get_initials(X,data,r,m)
  ini_alpha <- initials$alpha
  ini_beta <- initials$beta
  ini_delta <- initials$delta
  
  # size of parameters
  nalpha <- length(ini_alpha)
  nbeta <- length(ini_beta)
  ndelta <- length(ini_delta)

  # iteration counter
  iter <- 0
  # collect all costs
  costs <- c()
  # M-step update rate
  w <- 0.6 
  # indicator convergence
  isconverge <- FALSE; flag <- FALSE

  # EM algorithm loop
  while(!isconverge){

    # E-step
    pp <- p.matrix(X,ini_alpha,ini_beta,ini_delta,r,m) # risk probability matrix of two dieases
    tau.Y <- tau.Y1Y2(data,pp,Se,Sp,cj)$tau # conditional expectation
    costs <- c(costs,-E.loglik(ini_alpha,ini_delta,ini_beta,X,data,r,m,tau.Y)) # record history costs
    # M-step
    while(TRUE){
      
      # beta optimize by convolutional gradient 
      tem_beta <- optim(par=ini_beta, fn=l.beta, method='CG',
                        tem.alpha=ini_alpha,tem.delta=ini_delta,
                        X=X,data=data,r=r,m=m,tau.Y=tau.Y)$par
      # update beta with weight multiply difference
      tem_beta = ini_beta + w * (tem_beta-ini_beta)
      # w = 1.2; tem_beta = (1+w)*tem_beta - w*ini_beta

      # delta optimize
      tem_delta <- nlminb(start=ini_delta, objective=l.delta, 
                         lower = 0.001,upper=0.999,
                         tem.alpha=ini_alpha,ini.beta=tem_beta,tau.Y=tau.Y, 
                         X=X,data=data,r=r,m=m)$par
      
      # check convergence
      ini_cost <- -E.loglik(ini_alpha,ini_delta,ini_beta,X,data,r,m,tau.Y)
      tem_cost <- -E.loglik(ini_alpha,tem_delta,tem_beta,X,data,r,m,tau.Y)
      err <- abs(tem_cost-ini_cost)

      # iteration progress
      cat(paste0('estimated beta as: ')); cat(round(tem_beta,4),'\n')
      cat(paste0('estimated delta as: ')); cat(round(tem_delta,4),'\n')
      cat(paste0('consecutive costs difference as: ')); cat(round(err,4),'\n')
      flush.console()

      # reset the parameters
      ini_beta = tem_beta
      ini_delta = tem_delta

      # stopping rule
      if(flag | (err>1)){
        break
      }else{
        if(err<1e-2){
          flag=TRUE;
          break
        }
      }
    }

    # log transform to un-constrained
    ini_alpha <- nn_trans(ini_alpha, log_trans=TRUE)

    # spline coefficient alpha optimization
    tem_alpha <- optim(par=ini_alpha, fn=l.alpha,
                        ini.beta=ini_beta,ini.delta=ini_delta,
                        X=X,data=data,r=r,m=m,tau.Y=tau.Y)$par

    # exponential transform to non-negative
    ini_alpha <- nn_trans(ini_alpha)
    tem_alpha <- nn_trans(tem_alpha)

    # show iteration plots
    if(verbose){
      par(mfrow=c(1,2))
      u<-sort(X%*%c(1,tem_beta))
      alpha1 <- tem_alpha[1:(nalpha/2)]
      alpha2 <- tem_alpha[(nalpha/2+1):nalpha]
      plot(u,g(u,alpha1,r,m),type='l',col='gray',lty = 2,ylab=expression(hat(g)[1](u)))
      lines(u,links[[1]](u),type='l',col='red')
      plot(u,g(u,alpha2,r,m),type='l',col='gray',lty = 2,ylab=expression(hat(g)[2](u)))
      lines(u,links[[2]](u),type='l',col='red')
    }

    # update the loop
    iter <- iter + 1 
    ini_alpha <- tem_alpha
    ini_beta <- tem_beta
    ini_delta <- tem_delta

    # convergence status
    if(flag) isconverge <- TRUE
  }

  # return the estimates
  return(list(iterations=iter, costs=costs,
              alpha1=as.vector(tem_alpha)[1:(nalpha/2)],
              alpha2=as.vector(tem_alpha)[(nalpha/2+1):nalpha],
              beta=tem_beta,
              delta=tem_delta))
}