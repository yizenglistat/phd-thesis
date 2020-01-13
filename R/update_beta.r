fn_beta<- function(par,alpha,delta,X,data,ord,niknots,tau.Y){
  return(-E.loglik(alpha,delta,par,X,data,ord,niknots,tau.Y))
}

update_beta <- function(curr_alpha,curr_beta,curr_delta,w,optimizer,X,data,ord,niknots,tau.Y){
	# beta optimize by convolutional gradient 
	if(optimizer == "optim"){
		next_beta <- optim(par=curr_beta, fn=fn_beta, method='CG',
	                alpha=curr_alpha,delta=curr_delta,
	                X=X,data=data,ord=ord,niknots=niknots,tau.Y=tau.Y)$par
	}else if(optimizer == "nlminb"){
		next_beta <- nlminb(start=curr_beta, objective=fn_beta,
	                alpha=curr_alpha,delta=curr_delta,
	                X=X,data=data,ord=ord,niknots=niknots,tau.Y=tau.Y)$par
	}else{
		stop("optimizer not recoginized in update_beta function!")
	}
	
	# update beta with weight multiply difference
	next_beta <- curr_beta + w * (next_beta-curr_beta)
	#next_beta <- next_beta + w * (next_beta-curr_beta)

	return(next_beta)
}