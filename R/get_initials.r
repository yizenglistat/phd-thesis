"
This function generates initials for EM algorithm

Args:

  X       : features matrix
  data    : observed data matrix
  ord     : order of spline
  niknots : number of interior knots

return:
  
  A list of initials containing alpha, beta and delta
"
get_initials <- function(X, data, ord, niknots){
  # initial all zero for beta
  ini_beta <- rep(0,ncol(X)-1)
  # initial logisitc function as splines estimated functions (for both)
  u = seq(-8,8,length.out = 1000)
  ini_alpha = pnnls(Bsp(u,ord,niknots),u,1)$x; ini_alpha=c(ini_alpha,ini_alpha)
  # uniform initial delta, set edge to 0.001 for algorithm stability
  ini_delta <- runif(1,0.01,0.99)
  
  return(list(alpha=ini_alpha,beta=ini_beta,delta=ini_delta))
}


