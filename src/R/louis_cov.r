# Matrix of P(TEST|TRUEs), Six Cases 
# Output: 6 matrices, of which J's cjxcj matrices on the diagonal
TestCondMat<- function(p.hat,data,qY,qZ,N,cj,Se,Sp){
  # Special case when cj=2
  Ind<- out.prod(rep(1,N),rep(1,N),cj,1)
  p11.s0<- (1-Se[1])*Sp[2]*Ind 
  p22.s0<- Sp[1]*(1-Se[2])*Ind
  p11.s1<- qZ[,2]*out.prod(qY[,2],qY[,2],cj,1)
  p22.s1<- qZ[,3]*out.prod(qY[,3],qY[,3],cj,1)
  p33.s1<- qZ[,4]*out.prod(qY[,4],qY[,4],cj,1)
  p12.s1<- qZ[,4]*out.prod(qY[,2],qY[,3],cj,1)
  p13.s1<- qZ[,4]*out.prod(qY[,2],qY[,4],cj,1)
  p23.s1<- qZ[,4]*out.prod(qY[,3],qY[,4],cj,1)
  # Z=0
  Ind<- out.prod(rep(1,N),rep(1,N),cj,0)
  p11.0<- ((1-Se[1])*(1-Se[2])*Ind+(1-Se[1])*(Se[2]+Sp[2]-1)*prod.del2(p.hat[,1]+p.hat[,2],cj))+p11.s0
  p22.0<- ((1-Se[1])*(1-Se[2])*Ind+(1-Se[2])*(Se[1]+Sp[1]-1)*prod.del2(p.hat[,1]+p.hat[,3],cj))+p22.s0
  p33.0<- p12.0<- p13.0<- p23.0<- (1-Se[1])*(1-Se[2])*out.prod(rep(1,N),rep(1,N),cj,3) 
  # Z=1
  tem1 = prod.del2(qY[,1]*p.hat[,1]+qY[,2]*p.hat[,2],cj)
  tem2 = prod.del2(qY[,3]*p.hat[,3]+qY[,4]*p.hat[,4],cj)
  p11.1<- out.prod(qY[,2],qY[,2],cj,0)*
    (qZ[,2]*tem1+qZ[,4]*tem2)+p11.s1
  p22.1<- out.prod(qY[,3],qY[,3],cj,0)*(qZ[,3]*prod.del2(qY[,1]*p.hat[,1]+qY[,3]*p.hat[,3],cj)+qZ[,4]*prod.del2(qY[,2]*p.hat[,2]+qY[,4]*p.hat[,4],cj))+p22.s1
  temp = (qZ[,4]*prod.del2(qY[,1]*p.hat[,1]+qY[,2]*p.hat[,2]+qY[,3]*p.hat[,3]+qY[,4]*p.hat[,4],cj))
  p33.1<- out.prod(qY[,4],qY[,4],cj,0)*temp+p33.s1
  p12.1<- out.prod(qY[,2],qY[,3],cj,0)*temp+p12.s1
  p13.1<- out.prod(qY[,2],qY[,4],cj,0)*temp+p13.s1
  p23.1<- out.prod(qY[,3],qY[,4],cj,0)*temp+p23.s1
  # Merge
  Z = apply(data[,c(3,4)],1,max)
  p11<- p11.1*Z+p11.0*(1-Z)
  p22<- p22.1*Z+p22.0*(1-Z)
  p33<- p33.1*Z+p33.0*(1-Z)
  p12<- p12.1*Z+p12.0*(1-Z)
  p13<- p13.1*Z+p13.0*(1-Z)
  p23<- p23.1*Z+p23.0*(1-Z)
  return(list(p11=p11,p22=p22,p33=p33,p12=p12,p13=p13,p23=p23))
  return(list(p11=p11,p22=p22,p33=p33,p12=p12,p13=p13,p23=p23))
}

# Louis Method to compute the Co\textrm{Var}iance Matrix
LouisM <- function (data,X1,X2,start.beta1,start.beta2,start.delta,start.Se,start.Sp,N,cj){
  # Individual testing
  if(all(cj==1)){
    pp<- p.matrix(X1,X2,start.beta1,start.beta2,start.delta)
    eta.Y<-  Eta.Y1Y2(data,pp,start.Se,start.Sp,N,cj)$eta
    p1 = ncol(X1)-1; p2=ncol(X2)-1
    ppd<- pp.d(X1,X2,start.beta1,start.beta2,start.delta,p1+p2+3)
    ppd2<- pp.dd(X1,X2,start.beta1,start.beta2,start.delta,p1+p2+3)
    # calculate the expectation term
    l.d2<- 
      -t(sqrt(eta.Y[,1])/(pp[,1])*ppd$p00.d)%*%(sqrt(eta.Y[,1])/(pp[,1])*ppd$p00.d)-t(sqrt(eta.Y[,2])/(pp[,2])*ppd$p10.d)%*%(sqrt(eta.Y[,2])/(pp[,2])*ppd$p10.d)-
      t(sqrt(eta.Y[,3])/(pp[,3])*ppd$p01.d)%*%(sqrt(eta.Y[,3])/(pp[,3])*ppd$p01.d)-t(sqrt(eta.Y[,4])/(pp[,4])*ppd$p11.d)%*%(sqrt(eta.Y[,4])/(pp[,4])*ppd$p11.d)+
      Reduce("+",listOps(eta.Y[,1]/pp[,1],ppd2$p00.dd,"*"))+Reduce("+",listOps(eta.Y[,2]/pp[,2],ppd2$p10.dd,"*"))+
      Reduce("+",listOps(eta.Y[,3]/pp[,3],ppd2$p01.dd,"*"))+Reduce("+",listOps(eta.Y[,4]/pp[,4],ppd2$p11.dd,"*"))
    # calculate the "cov" term
    Q<- rbind(cbind(ppd$p10.d/(pp[,2])-ppd$p00.d/(pp[,1])), 
              cbind(ppd$p01.d/(pp[,3])-ppd$p00.d/(pp[,1])),
              cbind(ppd$p11.d/(pp[,4])-ppd$p00.d/(pp[,1])))
    # calculate the upper left 3 by 3 sub-block matrix of cov(V) 
    V11<- Diagonal(x=(eta.Y[,2]-eta.Y[,2]^2))
    V22<- Diagonal(x=(eta.Y[,3]-eta.Y[,3]^2))
    V33<- Diagonal(x=(eta.Y[,4]-eta.Y[,4]^2))
    V12<- Diagonal(x=(-eta.Y[,2]*eta.Y[,3]))
    V13<- Diagonal(x=(-eta.Y[,2]*eta.Y[,4]))
    V23<- Diagonal(x=(-eta.Y[,3]*eta.Y[,4]))
    COV<- t(Q)%*%rbind(cbind(V11,V12,V13),cbind(t(V12),V22,V23),cbind(t(V13),t(V23),V33))%*%Q
    # the hessian matrix
    H=l.d2+COV
  } else{
    p1 = ncol(X1)-1; p2 = ncol(X2)-1; J = length(cj); cj.vec = rep(cj,times=cj); c1 = 1*(cj.vec==1); Z = apply(data[,c(3,4)],1,max)
    pp<- p.matrix(X1,X2,start.beta1,start.beta2,start.delta)
    eta.Y<-  Eta.Y1Y2(data,pp,start.Se,start.Sp,N,cj)$eta
    eta.Z = eta_Z(data,pp,start.Se,start.Sp,N,cj)$eta
    pr.D<- Eta.Y1Y2(data,pp,start.Se,start.Sp,N,cj)$probs
    ppd<- pp.d(X1,X2,start.beta1,start.beta2,start.delta,p1+p2+3)
    ppd2<- pp.dd(X1,X2,start.beta1,start.beta2,start.delta,p1+p2+3)
    # calculate the expectation term
    h1.d2<- -t(sqrt(eta.Y[,1])/(pp[,1])*ppd$p00.d)%*%(sqrt(eta.Y[,1])/(pp[,1])*ppd$p00.d)-t(sqrt(eta.Y[,2])/(pp[,2])*ppd$p10.d)%*%(sqrt(eta.Y[,2])/(pp[,2])*ppd$p10.d)-
      t(sqrt(eta.Y[,3])/(pp[,3])*ppd$p01.d)%*%(sqrt(eta.Y[,3])/(pp[,3])*ppd$p01.d)-t(sqrt(eta.Y[,4])/(pp[,4])*ppd$p11.d)%*%(sqrt(eta.Y[,4])/(pp[,4])*ppd$p11.d)+
      Reduce("+",listOps(eta.Y[,1]/pp[,1],ppd2$p00.dd,"*"))+Reduce("+",listOps(eta.Y[,2]/pp[,2],ppd2$p10.dd,"*"))+
      Reduce("+",listOps(eta.Y[,3]/pp[,3],ppd2$p01.dd,"*"))+Reduce("+",listOps(eta.Y[,4]/pp[,4],ppd2$p11.dd,"*"))
    h2.d2<- diag(c(sum((1-c1)*(eta.Z[,2]+eta.Z[,4])/cj.vec*(data[,3]*(2*start.Se[1]-1)-start.Se[1]^2)+
                         (apply(cbind(Z,c1),1,max))*(eta.Y[,2]+eta.Y[,4])*(data[,1]*(2*start.Se[1]-1)-start.Se[1]^2))/(start.Se[1]^2*(1-start.Se[1])^2),
                   sum((1-c1)*(eta.Z[,3]+eta.Z[,4])/cj.vec*(data[,4]*(2*start.Se[2]-1)-start.Se[2]^2)+
                         (apply(cbind(Z,c1),1,max))*(eta.Y[,3]+eta.Y[,4])*(data[,2]*(2*start.Se[2]-1)-start.Se[2]^2))/(start.Se[2]^2*(1-start.Se[2])^2),
                   sum((1-c1)*(eta.Z[,1]+eta.Z[,3])/cj.vec*((1-data[,3])*(2*start.Sp[1]-1)-start.Sp[1]^2)+
                         (apply(cbind(Z,c1),1,max))*(eta.Y[,1]+eta.Y[,3])*((1-data[,1])*(2*start.Sp[1]-1)-start.Sp[1]^2))/(start.Sp[1]^2*(1-start.Sp[1])^2),
                   sum((1-c1)*(eta.Z[,1]+eta.Z[,2])/cj.vec*((1-data[,4])*(2*start.Sp[2]-1)-start.Sp[2]^2)+
                         (apply(cbind(Z,c1),1,max))*(eta.Y[,1]+eta.Y[,2])*((1-data[,2])*(2*start.Sp[2]-1)-start.Sp[2]^2))/(start.Sp[2]^2*(1-start.Sp[2])^2)))
    l.d2 = bdiag(h1.d2,h2.d2)    
    # calculate the "cov" term
    Q<- rbind(cbind2(ppd$p10.d/(pp[,2])-ppd$p00.d/(pp[,1]),sparseMatrix(dims=c(N,4),i={},j={})), 
              cbind2(ppd$p01.d/(pp[,3])-ppd$p00.d/(pp[,1]),sparseMatrix(dims=c(N,4),i={},j={})),
              cbind2(ppd$p11.d/(pp[,4])-ppd$p00.d/(pp[,1]),sparseMatrix(dims=c(N,4),i={},j={})), 
              1/(start.Se[1]*(1-start.Se[1]))*sparseMatrix(dims = c(N,p1+p2+7), i={1:N},j={rep(p1+p2+4,N)},x=rep(1,N)),
              1/(start.Se[2]*(1-start.Se[2]))*sparseMatrix(dims = c(N,p1+p2+7), i={1:N},j={rep(p1+p2+5,N)},x=rep(1,N)),
              1/(start.Sp[1]*(1-start.Sp[1]))*sparseMatrix(dims = c(N,p1+p2+7), i={1:N},j={rep(p1+p2+6,N)},x=rep(1,N)),
              1/(start.Sp[2]*(1-start.Sp[2]))*sparseMatrix(dims = c(N,p1+p2+7), i={1:N},j={rep(p1+p2+7,N)},x=rep(1,N)))
    # conditional prob of P(TESTs|TRUEs) for each individual
    qY = JointProb(data[,c(1,2)],start.Se,start.Sp)
    # conditional prob of P(TESTs|TRUEs) for each group
    qZ = JointProb(data[,c(3,4)],start.Se,start.Sp)
    # calculate the upper left 3 by 3 sub-block matrix of Cov(V) 
    pcM <- TestCondMat(pp,data,qY,qZ,N,cj,start.Se,start.Sp)
    V11 <- (pcM$p11/pr.D)*out.prod(pp[,2],pp[,2],cj,3)
    V22 <- (pcM$p22/pr.D)*out.prod(pp[,3],pp[,3],cj,3)
    V33 <- (pcM$p33/pr.D)*out.prod(pp[,4],pp[,4],cj,3)
    V12 <- (pcM$p12/pr.D)*out.prod(pp[,2],pp[,3],cj,3)
    V13 <- (pcM$p13/pr.D)*out.prod(pp[,2],pp[,4],cj,3)
    V23 <- (pcM$p23/pr.D)*out.prod(pp[,3],pp[,4],cj,3)
    diag(V11)<- eta.Y[,2]
    diag(V22)<- eta.Y[,3]
    diag(V33)<- eta.Y[,4]
    diag(V12)<- diag(V13)<- diag(V23)<- 0
    V11 <- V11-out.prod(eta.Y[,2],eta.Y[,2],cj,3)
    V22 <- V22-out.prod(eta.Y[,3],eta.Y[,3],cj,3)
    V33 <- V33-out.prod(eta.Y[,4],eta.Y[,4],cj,3)
    V12 <- V12-out.prod(eta.Y[,2],eta.Y[,3],cj,3)
    V13 <- V13-out.prod(eta.Y[,2],eta.Y[,4],cj,3)
    V23 <- V23-out.prod(eta.Y[,3],eta.Y[,4],cj,3)
    # calculate the rest of Cov(V) 
    Ind<- bdiag(lapply(cj, function(x){matrix(1, nrow=x, ncol=x)})) # index matrix 
    gama <- Gama(pp,cj)
    b <- GroupProb.del1(pp,data[,c(1,2)],start.Se,start.Sp,cj)
    V14<- as.vector((1-c1)*(data[,3]-start.Se[1])/cj.vec)*eta.Y[,2]*(eta.Z[,1]+eta.Z[,3])*Ind + as.vector((apply(cbind(Z,c1),1,max))*(data[,1]-start.Se[1]))*(V11+V13)
    tem01_1 = (as.vector(((1-start.Se[1])*(1-start.Se[2])*(gama[,2]+gama[,4])*pp[,3]*(1-Z)+(b[,2]+b[,4])*qZ[,4]*qY[,3]*pp[,3]*Z)/pr.D - eta.Y[,3]*(eta.Z[,2]+eta.Z[,4]))*Ind)
    V24<- as.vector((1-c1)*(data[,3]-start.Se[1])/cj.vec)*tem01_1 + as.vector((apply(cbind(Z,c1),1,max))*(data[,1]-start.Se[1]))*(V23+t(V12))
    tem11_1 = eta.Y[,4]*(eta.Z[,1]+eta.Z[,3])*Ind
    V34<- as.vector((1-c1)*(data[,3]-start.Se[1])/cj.vec)*tem11_1 + as.vector((apply(cbind(Z,c1),1,max))*(data[,1]-start.Se[1]))*(V33+t(V13))
    V44<- out.prod((1-c1)*(data[,3]-start.Se[1])/cj.vec*(eta.Z[,2]+eta.Z[,4]),(data[,3]-start.Se[1])/cj.vec*(eta.Z[,1]+eta.Z[,3]),cj,3) + 
      out.prod(apply(cbind(Z,c1),1,max)*(data[,1]-start.Se[1]),(data[,1]-start.Se[1]),cj,3)*(V11+V33+V13+t(V13))+
      out.prod((1-c1)*(data[,3]-start.Se[1])/cj.vec*(eta.Z[,1]+eta.Z[,3]),apply(cbind(Z,c1),1,max)*(data[,1]-start.Se[1])*(eta.Y[,2]+eta.Y[,4]),cj,3)+
      out.prod(apply(cbind(Z,c1),1,max)*(data[,1]-start.Se[1])*(eta.Y[,2]+eta.Y[,4]),(1-c1)*(data[,3]-start.Se[1])/cj.vec*(eta.Z[,1]+eta.Z[,3]),cj,3)
    tem10_2 = (as.vector(((1-start.Se[1])*(1-start.Se[2])*(gama[,3]+gama[,4])*pp[,2]*(1-Z)+(b[,3]+b[,4])*qZ[,4]*qY[,2]*pp[,2]*Z)/pr.D - eta.Y[,2]*(eta.Z[,3]+eta.Z[,4]))*Ind)
    V15<- as.vector((1-c1)*(data[,4]-start.Se[2])/cj.vec)*tem10_2+ as.vector((apply(cbind(Z,c1),1,max))*(data[,2]-start.Se[2]))*(V13+V12)
    V25<- as.vector((1-c1)*(data[,4]-start.Se[2])/cj.vec)*eta.Y[,3]*(eta.Z[,1]+eta.Z[,2])*Ind + as.vector((apply(cbind(Z,c1),1,max))*(data[,2]-start.Se[2]))*(V22+V23)
    tem11_2 = eta.Y[,4]*(eta.Z[,1]+eta.Z[,2])*Ind
    V35<- as.vector((1-c1)*(data[,4]-start.Se[2])/cj.vec)*tem11_2 + as.vector((apply(cbind(Z,c1),1,max))*(data[,2]-start.Se[2]))*(V33+t(V23))
    V45 = (eta.Z[,4]-(eta.Z[,2]+eta.Z[,4])*(eta.Z[,3]+eta.Z[,4]))*out.prod((1-c1)*(data[,3]-start.Se[1])/cj.vec,(data[,4]-start.Se[2])/cj.vec,cj,3) +
      out.prod(apply(cbind(Z,c1),1,max)*(data[,1]-start.Se[1]),(data[,2]-start.Se[2]),cj,3)*(V12+V13+t(V23)+V33)+
      out.prod((1-c1)*(data[,3]-start.Se[1])/cj.vec,apply(cbind(Z,c1),1,max)*(data[,2]-start.Se[2]),cj,3)*(t(tem01_1+tem11_1))+
      out.prod(apply(cbind(Z,c1),1,max)*(data[,1]-start.Se[1]),(1-c1)*(data[,4]-start.Se[2])/cj.vec,cj,3)*(tem10_2+tem11_2)
    V55<- out.prod((1-c1)*(data[,4]-start.Se[2])/cj.vec*(eta.Z[,3]+eta.Z[,4]),(data[,4]-start.Se[2])/cj.vec*(eta.Z[,1]+eta.Z[,2]),cj,3) + 
      out.prod((apply(cbind(Z,c1),1,max))*(data[,2]-start.Se[2]),(data[,2]-start.Se[2]),cj,3)*(V22+V33+V23+t(V23))+
      out.prod((1-c1)*(data[,4]-start.Se[2])/cj.vec*(eta.Z[,1]+eta.Z[,2]),apply(cbind(Z,c1),1,max)*(data[,2]-start.Se[2])*(eta.Y[,3]+eta.Y[,4]),cj,3)+
      out.prod(apply(cbind(Z,c1),1,max)*(data[,2]-start.Se[2])*(eta.Y[,3]+eta.Y[,4]),(1-c1)*(data[,4]-start.Se[2])/cj.vec*(eta.Z[,1]+eta.Z[,2]),cj,3)
    V16<- -(as.vector((1-c1)*(1-data[,3]-start.Sp[1])/cj.vec)*eta.Y[,2]*(eta.Z[,1]+eta.Z[,3])*Ind + as.vector((apply(cbind(Z,c1),1,max))*(1-data[,1]-start.Sp[1]))*(V11+V13))
    V26<- -(as.vector((1-c1)*(1-data[,3]-start.Sp[1])/cj.vec)*
              (as.vector(((1-start.Se[1])*(1-start.Se[2])*(gama[,2]+gama[,4])*pp[,3]*(1-Z)+(b[,2]+b[,4])*qZ[,4]*qY[,3]*pp[,3]*Z)/pr.D - eta.Y[,3]*(eta.Z[,2]+eta.Z[,4]))*Ind)+ 
              as.vector((apply(cbind(Z,c1),1,max))*(1-data[,1]-start.Sp[1]))*(V23+t(V12)))
    V36<- -(as.vector((1-c1)*(1-data[,3]-start.Sp[1])/cj.vec)*eta.Y[,4]*(eta.Z[,1]+eta.Z[,3])*Ind + as.vector((apply(cbind(Z,c1),1,max))*(1-data[,1]-start.Sp[1]))*(V33+t(V13)))
    V46<- -(out.prod((1-c1)*(data[,3]-start.Se[1])/cj.vec*(eta.Z[,2]+eta.Z[,4]),(1-data[,3]-start.Sp[1])/cj.vec*(eta.Z[,1]+eta.Z[,3]),cj,3) + 
              out.prod((apply(cbind(Z,c1),1,max))*(data[,1]-start.Se[1]),(1-data[,1]-start.Sp[1]),cj,3)*(V11+V33+V13+t(V13))+
              out.prod((1-c1)*(data[,3]-start.Se[1])/cj.vec*(eta.Z[,1]+eta.Z[,3]),apply(cbind(Z,c1),1,max)*(1-data[,1]-start.Sp[1])*(eta.Y[,2]+eta.Y[,4]),cj,3)+
              out.prod(apply(cbind(Z,c1),1,max)*(data[,1]-start.Se[1])*(eta.Y[,2]+eta.Y[,4]),(1-c1)*(1-data[,3]-start.Sp[1])/cj.vec*(eta.Z[,1]+eta.Z[,3]),cj,3))
    V56 = - ((eta.Z[,4]-(eta.Z[,3]+eta.Z[,4])*(eta.Z[,2]+eta.Z[,4]))*out.prod((1-c1)*(data[,4]-start.Se[2])/cj.vec,(1-data[,3]-start.Sp[1])/cj.vec,cj,3)+
               out.prod((1-c1)*(data[,4]-start.Se[2])/cj.vec,apply(cbind(Z,c1),1,max)*(1-data[,1]-start.Sp[1]),cj,3)*t(tem10_2+tem11_2)+
               out.prod(apply(cbind(Z,c1),1,max)*(data[,2]-start.Se[2]),(1-c1)*(1-data[,3]-start.Sp[1])/cj.vec,cj,3)*(tem01_1+tem11_1)+
               out.prod(apply(cbind(Z,c1),1,max)*(data[,2]-start.Se[2]),(1-data[,1]-start.Sp[1]),cj,3)*(t(V12)+t(V13)+V23+V33))
    V66<- out.prod((1-c1)*(1-data[,3]-start.Sp[1])/cj.vec*(eta.Z[,2]+eta.Z[,4]),(1-data[,3]-start.Sp[1])/cj.vec*(eta.Z[,1]+eta.Z[,3]),cj,3) + 
      out.prod((apply(cbind(Z,c1),1,max))*(1-data[,1]-start.Sp[1]),(1-data[,1]-start.Sp[1]),cj,3)*(V11+V33+V13+t(V13))+
      out.prod((1-c1)*(1-data[,3]-start.Sp[1])/cj.vec*(eta.Z[,1]+eta.Z[,3]),apply(cbind(Z,c1),1,max)*(1-data[,1]-start.Sp[1])*(eta.Y[,2]+eta.Y[,4]),cj,3)+
      out.prod(apply(cbind(Z,c1),1,max)*(1-data[,1]-start.Sp[1])*(eta.Y[,2]+eta.Y[,4]),(1-c1)*(1-data[,3]-start.Sp[1])/cj.vec*(eta.Z[,1]+eta.Z[,3]),cj,3)
    V17<- -(as.vector((1-c1)*(1-data[,4]-start.Sp[2])/cj.vec)*
              (as.vector(((1-start.Se[1])*(1-start.Se[2])*(gama[,3]+gama[,4])*pp[,2]*(1-Z)+(b[,3]+b[,4])*qZ[,4]*qY[,2]*pp[,2]*Z)/pr.D - eta.Y[,2]*(eta.Z[,3]+eta.Z[,4]))*Ind)+
              as.vector((apply(cbind(Z,c1),1,max))*(1-data[,2]-start.Sp[2]))*(V13+V12))
    V27<- -(as.vector((1-c1)*(1-data[,4]-start.Sp[2])/cj.vec)*eta.Y[,3]*(eta.Z[,1]+eta.Z[,2])*Ind + as.vector((apply(cbind(Z,c1),1,max))*(1-data[,2]-start.Sp[2]))*(V22+V23))
    V37<- -(as.vector((1-c1)*(1-data[,4]-start.Sp[2])/cj.vec)*eta.Y[,4]*(eta.Z[,1]+eta.Z[,2])*Ind + as.vector((apply(cbind(Z,c1),1,max))*(1-data[,2]-start.Sp[2]))*(V33+t(V23)))
    V47 = - ((eta.Z[,4]-(eta.Z[,2]+eta.Z[,4])*(eta.Z[,3]+eta.Z[,4]))*out.prod((1-c1)*(data[,3]-start.Se[1])/cj.vec,(1-data[,4]-start.Sp[2])/cj.vec,cj,3) +
               out.prod((1-c1)*(data[,3]-start.Se[1])/cj.vec,apply(cbind(Z,c1),1,max)*(1-data[,2]-start.Sp[2]),cj,3)*t(tem01_1+tem11_1)+
               out.prod(apply(cbind(Z,c1),1,max)*(data[,1]-start.Se[1]),(1-c1)*(1-data[,4]-start.Sp[2])/cj.vec,cj,3)*(tem10_2+tem11_2)+
               out.prod(apply(cbind(Z,c1),1,max)*(data[,1]-start.Se[1]),(1-data[,2]-start.Sp[2]),cj,3)*(V12+V13+t(V23)+V33))
    V57<- -(out.prod((1-c1)*(data[,4]-start.Se[2])/cj.vec*(eta.Z[,3]+eta.Z[,4]),(1-data[,4]-start.Sp[2])*(eta.Z[,1]+eta.Z[,2])/cj.vec,cj,3) +
              out.prod((apply(cbind(Z,c1),1,max))*(data[,2]-start.Se[2]),(1-data[,2]-start.Sp[2]),cj,3)*(V22+V33+V23+t(V23))+
              out.prod((1-c1)*(data[,4]-start.Se[2])/cj.vec*(eta.Z[,1]+eta.Z[,2]),(apply(cbind(Z,c1),1,max))*(1-data[,2]-start.Sp[2])*(eta.Y[,3]+eta.Y[,4]),cj,3)+
              out.prod((apply(cbind(Z,c1),1,max))*(data[,2]-start.Se[2])*(eta.Y[,3]+eta.Y[,4]),(1-c1)*(1-data[,4]-start.Sp[2])/cj.vec*(eta.Z[,1]+eta.Z[,2]),cj,3))
    V67 = (eta.Z[,4]-(eta.Z[,2]+eta.Z[,4])*(eta.Z[,3]+eta.Z[,4]))*out.prod((1-c1)*(1-data[,3]-start.Sp[1])/cj.vec,(1-data[,4]-start.Sp[2])/cj.vec,cj,3) +
      out.prod((1-c1)*(1-data[,3]-start.Sp[1])/cj.vec,apply(cbind(Z,c1),1,max)*(1-data[,2]-start.Sp[2]),cj,3)*t(tem01_1+tem11_1)+
      out.prod(apply(cbind(Z,c1),1,max)*(1-data[,1]-start.Sp[1]),(1-c1)*(1-data[,4]-start.Sp[2])/cj.vec,cj,3)*(tem10_2+tem11_2)+
      out.prod(apply(cbind(Z,c1),1,max)*(1-data[,1]-start.Sp[1]),(1-data[,2]-start.Sp[2]),cj,3)*(V12+V13+t(V23)+V33)
    V77 = out.prod((1-c1)*(1-data[,4]-start.Sp[2])/cj.vec*(eta.Z[,3]+eta.Z[,4]),(1-data[,4]-start.Sp[2])/cj.vec*(eta.Z[,1]+eta.Z[,2]),cj,3) + 
      out.prod((apply(cbind(Z,c1),1,max))*(1-data[,2]-start.Sp[2]),(1-data[,2]-start.Sp[2]),cj,3)*(V22+V33+V23+t(V23))+
      out.prod((1-c1)*(1-data[,4]-start.Sp[2])/cj.vec*(eta.Z[,1]+eta.Z[,2]),apply(cbind(Z,c1),1,max)*(1-data[,2]-start.Sp[2])*(eta.Y[,3]+eta.Y[,4]),cj,3)+
      out.prod(apply(cbind(Z,c1),1,max)*(1-data[,2]-start.Sp[2])*(eta.Y[,3]+eta.Y[,4]),(1-c1)*(1-data[,4]-start.Sp[2])/cj.vec*(eta.Z[,1]+eta.Z[,2]),cj,3)
    
    COV<- t(Q)%*%rbind(cbind(V11,V12,V13,V14,V15,V16,V17),cbind(t(V12),V22,V23,V24,V25,V26,V27),
                       cbind(t(V13),t(V23),V33,V34,V35,V36,V37),cbind(t(V14),t(V24),t(V34),V44,V45,V46,V47),
                       cbind(t(V15),t(V25),t(V35),t(V45),V55,V56,V57),cbind(t(V16),t(V26),t(V36),t(V46),t(V56),V66,V67),
                       cbind(t(V17),t(V27),t(V37),t(V47),t(V57),t(V67),V77))%*%Q
    
    # Hessian matrix
    H=l.d2+COV
  }
  return(H)
}