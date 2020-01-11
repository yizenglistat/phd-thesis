#*******************************************************************************#
#               Supporting functions for convenience purpose
#*******************************************************************************#
nn_trans <- function(vec,log_trans=FALSE){
  nvec = length(vec)
  if(log_trans){
    vec[2:(nvec/2)]<-log(vec[2:(nvec/2)])
    vec[(nvec/2+2):nvec]<-log(vec[(nvec/2+2):nvec])
  }else{
    vec[2:(nvec/2)]<-exp(vec[2:(nvec/2)])
    vec[(nvec/2+2):nvec]<-exp(vec[(nvec/2+2):nvec])
  }
  return(vec)
}

# list operation 
listOps <- function(a,b,Ops){
  return(mapply(Ops,a,b,SIMPLIFY = FALSE))
}

# Split vector into groups
# Output: a list of length = J, each element represents one group
Split <- function(a,cj){
  ind1<- c(1,cumsum(cj)[-length(cj)]+1)
  ind2<- ind1+cj-1
  return(lapply(data.frame(rbind(ind1,ind2)), function(x){a[x[1]:x[2]]}))
}

# Vectors outer product within each group 
# Input: vector a, vector b, vector cj, int special
#   cj = group sizes, a vector of length = J
#   special = 1/0 indicates special case when cj is 2 or not. 
#   special = 3 means we don't consider the special case
# Output: a block diagonal matrix, each block is a cj by cj matrix, 
#         the outer product of a*b of each group.
out.prod<- function(a,b,cj,special){
  cj.vec<- rep(cj,times=cj) # length = N
  if(special==1){
    a=(cj.vec==2)*a
    b=(cj.vec==2)*b
  } else if(special==0){
    a=(cj.vec!=2)*a
    b=(cj.vec!=2)*b
  } 
  if(all(a==b)){
    op <- bdiag(lapply(Split(a,cj),function(x){x%*%t(x)}))
  }
  else{
    op <- bdiag(mapply(function(x,y){x%*%t(y)},
                       Split(a,cj),Split(b,cj),SIMPLIFY = FALSE))
  }
  return(op)
}

# Input: vector a, vector cj
# Output: a vector of length = N, each entry indicates the group product  
prod.gr <- function(a,cj){
  prod.g <- unlist(lapply(Split(a,cj),function(x){prod(x)}))
  prod.g <- rep(prod.g,times=cj)
  return(prod.g)
}

# Input: vector a, vector cj
# Output: a vector of length = N, each entry indcicates the group product without itself. 
prod.del1 <- function(a,cj){
  prod.del <- unlist(lapply(Split(a,cj),function(x){prod(x)/x}))
  return(prod.del)
}

# Input: vector a, vector cj
# Output: a block diagonal matrix, where each block is a cj by cj matrix, 
#         indicating group multiplication without two subjects.
prod.del2<- function(a,cj){
  prod.del <- bdiag(lapply(Split(a,cj),function(x){prod(x)/(x%*%t(x))}))
  return(prod.del)
}