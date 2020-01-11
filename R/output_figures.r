output_figures <- function(curr_alpha, curr_beta, ord, niknots, links){
	par(mfrow=c(1,2))
	u<-sort(X%*%c(1,curr_beta))
	nalpha <- length(curr_alpha)
	alpha1 <- curr_alpha[1:(nalpha/2)]
	alpha2 <- curr_alpha[(nalpha/2+1):nalpha]
	plot(u,g(u,alpha1,ord,niknots),type='l',col='gray',lty = 2, ylab=expression(hat(g)[1](u)))
	lines(u,links[[1]](u),type='l',col='red')
	plot(u,g(u,alpha2,ord,niknots),type='l',col='gray',lty = 2, ylab=expression(hat(g)[2](u)))
	lines(u,links[[2]](u),type='l',col='red')
}