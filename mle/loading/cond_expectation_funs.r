#***********************************************************************#
#                   2.5 Conditioanl Expectations functions
#***********************************************************************#
# prevalence matrix 
# Input: 
#   X = covariates information; 
#   beta, alpha, delta, p, m = parameters
# Output: N by 4 matrix, of which columns are p_00, p_10, p_01, p_11 
p.matrix<- function(X,ini.alpha,ini.beta,ini.delta,r,m){
  N<-nrow(X)
  pp<- matrix(NA,N,4)
  
  u<-X%*%c(1,ini.beta)
  #u<-X%*%beta
  
  pp[,4]<- gumbel(g(u,ini.alpha[1:(r+m)],r,m),g(u,ini.alpha[(r+m+1):(2*(r+m))],r,m),ini.delta)
  #pp[,4]<- g(u,ini.alpha[1:(r+m)],r,m)*g(u,ini.alpha[(r+m+1):(2*(r+m))],r,m)
  pp[,2]<- g(u,ini.alpha[1:(r+m)],r,m)-pp[,4]
  pp[,3]<- g(u,ini.alpha[(r+m+1):(2*(r+m))],r,m)-pp[,4]
  pp[,1]<- 1-pp[,2]-pp[,3]-pp[,4]
  
  pp[pp<=1e-6]=1e-6
  
  return(pp)
}

# Prob of TRUEs for the jth group without considering i,jth individual
Gama<- function(p.mat,cj){
  gama<- matrix(NA,sum(cj),4)
  gama[,1]<- prod.del1(p.mat[,1],cj)
  gama[,2]<- prod.del1((p.mat[,1]+p.mat[,2]),cj) - gama[,1]
  gama[,3]<- prod.del1((p.mat[,1]+p.mat[,3]),cj) - gama[,1]
  gama[,4]<- 1-gama[,1]-gama[,2]-gama[,3]
  return(gama)
}

# Prob of group TESTs given individual TRUEs
JointProb<- function(data,Se,Sp){
  q<- matrix(NA,dim(data)[1],4)
  q[,1]<- Sp[1]^(1-data[,1])*(1-Sp[1])^data[,1]*Sp[2]^(1-data[,2])*(1-Sp[2])^data[,2]
  q[,2]<- Se[1]^data[,1]*(1-Se[1])^(1-data[,1])*Sp[2]^(1-data[,2])*(1-Sp[2])^data[,2]
  q[,3]<- Sp[1]^(1-data[,1])*(1-Sp[1])^data[,1]*Se[2]^data[,2]*(1-Se[2])^(1-data[,2])
  q[,4]<- Se[1]^data[,1]*(1-Se[1])^(1-data[,1])*Se[2]^data[,2]*(1-Se[2])^(1-data[,2])
  return(q)
}

# Joint prob of group TESTs and individual TRUEs for the jth group without considering i,jth individual
GroupProb.del1<- function(p.mat,data,Se,Sp,cj){
  q<-JointProb(data,Se,Sp)
  b<-matrix(NA,sum(cj),4)
  b[,1]<- prod.del1(q[,1]*p.mat[,1],cj)
  b[,2]<- prod.del1(q[,1]*p.mat[,1]+q[,2]*p.mat[,2],cj)-b[,1]
  b[,3]<- prod.del1(q[,1]*p.mat[,1]+q[,3]*p.mat[,3],cj)-b[,1]
  b[,4]<- prod.del1(rowSums(q*p.mat),cj)-b[,1]-b[,2]-b[,3]
  return(b)
}

# conditional expectations with interactions of Y1 and Y2
# or conditional probability pr(tilde Yij1,tilde Yij2 | Yij1, Yij2)
# output is a N by 4 matrix with each row being (tau_{ij00},tau_{ij10},tau_{ij01},tau_{ij11})
tau.Y1Y2<- function(data,p.hat,Se,Sp,cj){
  # individual testing 
  if(all(cj==1)){
    qY<- JointProb(data[,c(1,2)],Se,Sp)
    w<- cbind(qY[,1]*p.hat[,1],qY[,2]*p.hat[,2],qY[,3]*p.hat[,3],qY[,4]*p.hat[,4])
    tau <- w/rowSums(w)
    probs = rowSums(w)
  } else{ # group testing
    ## special case cj=1
    J<- length(cj)
    cj.s<- rep(ifelse(cj==1,1,0),times=cj)   # length(cj.s) = N 
    qZ<- JointProb(data[,c(3,4)],Se,Sp)
    w.s<- cbind(qZ[,1]*p.hat[,1],qZ[,2]*p.hat[,2],qZ[,3]*p.hat[,3],qZ[,4]*p.hat[,4])
    tau.s <- w.s/rowSums(w.s)
    ### Zj=0
    gama<-Gama(p.hat,cj)
    w0<-matrix(NA,N,4)
    w0[,1]<- (Sp[1]*Sp[2]*gama[,1]+(1-Se[1])*Sp[2]*gama[,2]+Sp[1]*(1-Se[2])*gama[,3]+(1-Se[1])*(1-Se[2])*gama[,4])*p.hat[,1]
    w0[,2]<- (Sp[2]*(gama[,1]+gama[,2])+(1-Se[2])*(gama[,3]+gama[,4]))*(1-Se[1])*p.hat[,2]
    w0[,3]<- (Sp[1]*(gama[,1]+gama[,3])+(1-Se[1])*(gama[,2]+gama[,4]))*(1-Se[2])*p.hat[,3]
    w0[,4]<- (1-Se[1])*(1-Se[2])*p.hat[,4]
    tau0<-w0/rowSums(w0)
    ### Zj=1
    qY<- JointProb(data[,c(1,2)],Se,Sp)
    b<- GroupProb.del1(p.hat,data[,c(1,2)],Se,Sp,cj)
    w1<- matrix(NA,N,4)
    w1[,1]<- (b[,1]*qZ[,1]+b[,2]*qZ[,2]+b[,3]*qZ[,3]+b[,4]*qZ[,4])*qY[,1]*p.hat[,1]
    w1[,2]<- ((b[,1]+b[,2])*qZ[,2]+(b[,3]+b[,4])*qZ[,4])*qY[,2]*p.hat[,2]
    w1[,3]<- ((b[,1]+b[,3])*qZ[,3]+(b[,2]+b[,4])*qZ[,4])*qY[,3]*p.hat[,3]
    w1[,4]<- (b[,1]+b[,2]+b[,3]+b[,4])*qZ[,4]*qY[,4]*p.hat[,4]
    tau1<-w1/rowSums(w1)
    ### merge
    probs<- (1-cj.s)*(rowSums(w1)*(data[,3]+data[,4]-data[,3]*data[,4])+rowSums(w0)*(1-data[,3])*(1-data[,4]))+cj.s*rowSums(w.s)
    tau<- (1-cj.s)*(tau1*(data[,3]+data[,4]-data[,3]*data[,4])+tau0*(1-data[,3])*(1-data[,4]))+cj.s*tau.s
  }
  return(list(tau=tau,probs=probs))
}