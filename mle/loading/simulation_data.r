#***********************************************************************#
#             Generate true prevalence probability
#***********************************************************************#
# true probability
p.matrix.true<- function(X,beta,link1,link2,delta){
  # beta <- beta/as.numeric(sqrt(t(beta)%*%beta))
  beta <- c(1,beta)
  
  N<-nrow(X)
  pp<- matrix(NA,N,4)
  
  u<-X%*%beta
  #u <- scale(u,center = min(u),scale = max(u)-min(u))
  
  pp[,4]<- gumbel(link1(u),link2(u),delta)
  pp[,2]<- link1(u)-pp[,4]
  pp[,3]<- link2(u)-pp[,4]
  pp[,1]<- 1-pp[,2]-pp[,3]-pp[,4]
  return(pp)
}

#***********************************************************************#
#             Generate simulation data
#***********************************************************************#
SIM<- function(X,Y,N,cj,Se,Sp){
  # Generate true pervalence prob matrix

  # create groups
  #J = ceiling(N/c)
  #cj = c(rep(c-1,100),rep(c,J-1)	,N-(J-1)*c)					
  J <- length(cj)
  ind1<- c(1,cumsum(cj)[-J]+1)
  ind2<- ind1+cj-1
  # simulate the testing outcomes
  data<-matrix(NA,nrow=N,ncol=4)
  num.test<-0
  for(j in 1:J){
    ind<-ind1[j]:ind2[j]
    if(length(ind)==1){
      prob<- 1-Sp+Y[ind,]*(Se+Sp-1)
    }
    if(length(ind)>1){
      prob<- 1-Sp+apply(Y[ind,],2,max)*(Se+Sp-1)
    }
    Z1<-rbinom(1,1,prob[1])
    Z2<-rbinom(1,1,prob[2])
    data[ind,3]<-Z1
    data[ind,4]<-Z2
    num.test<-num.test+1
    if(Z1==0 & Z2==0){
      data[ind,c(1,2)]<-0
    }
    if(length(ind)==1){
      data[ind,c(1,2)]<-c(Z1,Z2)
    }
    if(length(ind)>1){
      if(Z1==1 | Z2==1){
        for(i in ind){
          prob<- 1-Sp+Y[i,]*(Se+Sp-1)
          Y1<-rbinom(1,1,prob[1])
          Y2<-rbinom(1,1,prob[2])
          num.test<-num.test+1
          data[i,c(1,2)]<-c(Y1,Y2)
        }
      }
    }
  }
  return(list(cj=cj, num.test=num.test,data=data, Y=Y))
}

SIM_true <- function(X,beta,links,delta){
  # generate true prevalence probability 
  p.true<-p.matrix.true(X,beta,links[[1]],links[[2]],delta)
  # Generate true true responses
  Y <- matrix(NA,nrow=N,ncol=2)
  for(i in 1:N){
    stat <- rmultinom(1,1,p.true[i,])
    if(stat[1]==1){Y[i,]<-c(0,0)}
    if(stat[2]==1){Y[i,]<-c(1,0)}
    if(stat[3]==1){Y[i,]<-c(0,1)}
    if(stat[4]==1){Y[i,]<-c(1,1)}
  }
  return(Y)
}