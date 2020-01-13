output_figure <- function(curr_alpha, curr_beta, ord, niknots, links){
	par(mfrow=c(2,2))
	u<-sort(X%*%c(1,curr_beta))
	nalpha <- length(curr_alpha)
	alpha1 <- curr_alpha[1:(nalpha/2)]
	alpha2 <- curr_alpha[(nalpha/2+1):nalpha]

	eta1 <- function(u) log(links[[1]](u)/(1-links[[1]](u)))
	eta2 <- function(u) log(links[[2]](u)/(1-links[[2]](u)))
	etas <- list(eta1,eta2)

	plot(u,g(u,alpha1,ord,niknots),type='l',col='gray',lty = 2, ylab=expression(hat(g)[1](u)))
	lines(u,links[[1]](u),type='l',col='green')
	plot(u,g(u,alpha2,ord,niknots),type='l',col='gray',lty = 2, ylab=expression(hat(g)[2](u)))
	lines(u,links[[2]](u),type='l',col='green')

	plot(u,eta.Bsp(u,alpha1,ord,niknots),type='l',col='gray',lty = 2, ylab=expression(hat(eta)[1](u)))
	lines(u,etas[[1]](u),type='l',col='green')
	plot(u,eta.Bsp(u,alpha2,ord,niknots),type='l',col='gray',lty = 2, ylab=expression(hat(eta)[2](u)))
	lines(u,etas[[2]](u),type='l',col='green')
}