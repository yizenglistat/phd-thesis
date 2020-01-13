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
  isfull  : TRUE or FALSE, overwrite former iterations or show all iterations in the console
  seed    : seed state number

Return:

  A list of estimates and costs iteration history

"
mle <- function(X, cj, data, Se, Sp, ord, niknots, verbose=TRUE, isfull=FALSE, seed){
  
  # initial setting  
  initials <- get_initials(X,data,ord,niknots)
  curr_alpha <- initials$alpha
  curr_beta <- initials$beta
  curr_delta <- initials$delta
  
  # size of parameters
  nalpha <- length(curr_alpha)
  nbeta <- length(curr_beta)
  ndelta <- length(curr_delta)

  # iteration counter
  iter <- 0
  # initial cost
  pp <- prob_mat(X,curr_alpha,curr_beta,curr_delta,ord,niknots) # risk probability matrix of two dieases
  tau.Y <- tau.Y1Y2(data,pp,Se,Sp,cj)$tau # conditional expectation
  costs <- -E.loglik(curr_alpha,curr_delta,curr_beta,X,data,ord,niknots,tau.Y)
  # M-step update rate
  w <- 0.75
  # indicator convergence
  isconverge <- FALSE; err <- NA
  # save EM output
  output_file = paste0('./output/simulation/output',seed,'.md')
  figure_file = paste0('./output/simulation/figures/figure',seed,'.png');
  if (file.exists(output_file)) file.remove(output_file) # delete file if it exists
  if (file.exists(figure_file)) file.remove(figure_file) # delete file if it exists
  sep_lines <- paste0(paste0(rep('-',20+nbeta*7+(nbeta-1)),collapse =''))
  header <- output_header(N, ord, niknots, Se, Sp, beta, delta)
  cat('```r',sep_lines,paste0(paste(rep('',(nchar(sep_lines)-7)/2+1),sep=' ',collapse=' '),'Setting'),
    header, paste0(paste(rep('',(nchar(sep_lines)-24)/2+1),sep=' ',collapse=' '),'Accelerated EM Algorithm'), 
    sep='\n',file=output_file,append=TRUE)
  href <- 'source at https://github.com/yizenglistat/regression-supervised-multiple-infection-group-testing'
  # EM algorithm loop
  while(!isconverge){

    # nested EM for regression parameters
    while(TRUE){
      
      pp <- prob_mat(X,curr_alpha,curr_beta,curr_delta,ord,niknots) # risk probability matrix of two dieases
      tau.Y <- tau.Y1Y2(data,pp,Se,Sp,cj)$tau # conditional expectation

      # update beta
      next_beta <- update_beta(curr_alpha=curr_alpha,curr_beta=curr_beta,curr_delta=curr_delta, w=w, optimizer="optim",
                                X=X,data=data,ord=ord,niknots=niknots,tau.Y=tau.Y)
      next_delta <- update_delta(curr_alpha=curr_alpha,curr_beta=next_beta,curr_delta=curr_delta,
                                X=X,data=data,ord=ord,niknots=niknots,tau.Y=tau.Y)
      # check convergence
      curr_cost <- -E.loglik(curr_alpha,curr_delta,curr_beta,X,data,ord,niknots,tau.Y)
      next_cost <- -E.loglik(curr_alpha,next_delta,next_beta,X,data,ord,niknots,tau.Y)
      nest_err  <- curr_cost-next_cost

      # update beta if needed
      if(nest_err<0){
        next_beta <- update_beta(curr_alpha=curr_alpha,curr_beta=curr_beta,curr_delta=curr_delta, w=w, optimizer="nlminb",
                                X=X,data=data,ord=ord,niknots=niknots,tau.Y=tau.Y)
        next_delta <- update_delta(curr_alpha=curr_alpha,curr_beta=next_beta,curr_delta=curr_delta,
                                X=X,data=data,ord=ord,niknots=niknots,tau.Y=tau.Y)
      }
      
      # check convergence
      curr_cost <- -E.loglik(curr_alpha,curr_delta,curr_beta,X,data,ord,niknots,tau.Y)
      next_cost <- -E.loglik(curr_alpha,next_delta,next_beta,X,data,ord,niknots,tau.Y)
      nest_err <- curr_cost-next_cost

      nest_beta_err <- sqrt(sum((next_beta-curr_beta)^2))

      # iteration progress
      body <- output_body(costs, err, nest_err, nest_beta_err, next_beta, next_delta)
      if(verbose){
        if(!isfull) cat("\014")
        cat(sep_lines, green('setting'),header,green('accelerated EM algorithm'),body,sep_lines,'more detail at README.md',sep='\n')
      }
      cat(body,sep='\n',file=output_file,append=TRUE)

      # update nested EM loop
      if(nest_err<0){
        p <- runif(1); 
        curr_beta = next_beta * (p<0) + curr_beta*(p>=0)
        curr_delta = next_delta * (p<0) + curr_delta*(p>=0)
      }else{
        curr_beta = next_beta
        curr_delta = next_delta
      }

      # stopping rule
      if( (nest_err>1e0) | (abs(nest_err)<1e-1) ){
        break
      }
    }

    # nested E-step
    pp <- prob_mat(X,curr_alpha,curr_beta,curr_delta,ord,niknots) # risk probability matrix of two dieases
    tau.Y <- tau.Y1Y2(data,pp,Se,Sp,cj)$tau # conditional expectation
    
    # nested M-step update alpha
    next_alpha <- update_alpha(curr_alpha=curr_alpha,curr_beta=curr_beta,curr_delta=curr_delta,
                                X=X,data=data,ord=ord,niknots=niknots,tau.Y=tau.Y)

    # update nested EM loop
    iter <- iter + 1 
    curr_alpha <- next_alpha

    # show iteration plots
    if(verbose){
      if(isconverge) png(filename=figure_file) # save plots when converged
      output_figure(curr_alpha, curr_beta, ord, niknots, links)
      if(isconverge) dev.off()
    }

    # record costs history
    costs <- c(costs,-E.loglik(curr_alpha,curr_delta,curr_beta,X,data,ord,niknots,tau.Y)) 
    if(length(costs)>1){
      err <- (tail(costs,2)[1]-tail(costs,1))/tail(costs,2)[1] # change rate (percent)
      if((err>0) & (err<1e-2)) isconverge <- TRUE
    }
  }

  # output
  cat("\014")
  tail <- output_tail(costs, err, nest_err, nest_beta_err, curr_beta, curr_delta)
  cat(tail,'Done!',sep='\n')
  cat(tail,'```',href,sep='\n',file=output_file,append=TRUE)

  # return the estimates
  return(list(iterations=iter, costs=costs,
              alpha1=as.vector(curr_alpha)[1:(nalpha/2)],
              alpha2=as.vector(curr_alpha)[(nalpha/2+1):nalpha],
              beta=curr_beta,
              delta=curr_delta))
}