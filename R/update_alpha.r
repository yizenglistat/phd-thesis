"

Args:

Return:
	updated alpha vector

"

fn_alpha<- function(par,beta,delta,X,data,ord,niknots,tau.Y){
  # transform to nonnegative parameters
  par <- nn_trans(par)
  return(-E.loglik(par,delta,beta,X,data,ord,niknots,tau.Y))
}

update_alpha <- function(curr_alpha,curr_beta,curr_delta,X,data,ord,niknots,tau.Y){

	# log transform to un-constrained
	curr_alpha <- nn_trans(curr_alpha, log_trans=TRUE) 
	# spline coefficient alpha optimization
	next_alpha <- optim(par=curr_alpha, fn=fn_alpha, 
	                    beta=curr_beta,delta=curr_delta,
	                    X=X,data=data,ord=ord,niknots=niknots,tau.Y=tau.Y)$par

	# exponential transform to non-negative
	curr_alpha <- nn_trans(curr_alpha)
	next_alpha <- nn_trans(next_alpha)

	return(next_alpha)
}
