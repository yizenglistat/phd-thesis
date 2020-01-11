#***********************************************************************#
#                       Generized link functions
#***********************************************************************#
# inverse of logit function
logistic <- function(u){return(1/(1+exp(-u)))}

# generilized link function, g_k in our case, depend on alpha_k
g <- function(u,alphak,r,m){return(logistic(eta.Bsp(u,alphak,r,m)))}

# Gumbel copula
gumbel <- function(u,v,delta){
  return(exp(-((-log(u))^(1/delta)+(-log(v))^(1/delta))^(delta)))
}
