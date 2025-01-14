fn_delta<- function(par,alpha,beta,X,data,ord,niknots,tau.Y){
  return(-E.loglik(alpha,par,beta,X,data,ord,niknots,tau.Y))
}

update_delta <- function(curr_alpha,curr_beta,curr_delta,w,X,data,ord,niknots,tau.Y){
	# delta optimize
	next_delta <- nlminb(start=curr_delta, objective=fn_delta, 
	                  lower=0.001,upper=0.999,
	                  alpha=curr_alpha,beta=curr_beta,tau.Y=tau.Y, 
	                  X=X,data=data,ord=ord,niknots=niknots)$par
	temp_delta <- curr_delta + w * (next_delta-curr_delta)
	if((temp_delta>=1) | (temp_delta<=0)) temp_delta <- next_delta;
	next_delta <- temp_delta
	return(next_delta)
}