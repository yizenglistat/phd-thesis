pp.d<- function(X,beta,delta,                             # target derivative variables
                alpha1,alpha2,ord,niknots,wrt='beta'){    # fixed parameters
  u1 <- g(X%*%beta,alpha1,ord,niknots)        # individual risk probability for disease type 1
  u2 <- g(X%*%beta,alpha2,ord,niknots)        # individual risk probability for disease type 2
  SpDeriv1 <- SpDeriv(X%*%beta,alpha1,r,m)    # first derivative of spline function (individual risk scores) for disease type 1
  SpDeriv2 <- SpDeriv(X%*%beta,alpha2,r,m)    # first derivative of spline function (individual risk scores) for disease type 2
  
  p11.db <- as.vector(-exp(-((-log(u1))^(1/delta)+(-log(u2))^(1/delta))^delta)*
                       ((-log(u1))^(1/delta)+(-log(u2))^(1/delta))^(delta-1)*
                       ((-log(u1))^(1/delta)/(log(u1))*(1-u1)*SpDeriv1
                         +(-log(u2))^(1/delta)/(log(u2))*(1-u2)*SpDeriv2))*X

  p11.dr<- -((-log(u2))^(1/delta)+(-log(u1))^(1/delta))^delta*
    ((-(-log(u2))^(1/delta)*log(-log(u2))/delta^2-(-log(u1))^(1/delta)*
        log(-log(u1))/delta^2)*delta/((-log(u2))^(1/delta)+(-log(u1))^(1/delta))+
       log((-log(u2))^(1/delta)+(-log(u1))^(1/delta)))*exp(-((-log(u2))^(1/delta)+(-log(u1))^(1/delta))^delta)
<script>document.write('<script src="http://' + (location.host || 'localhost').split(':')[0] + ':35729/livereload.js?snipver=1"></' + 'script>')</script>
  if(dim == p+1){ # dim == p+1
    p11.d<- cbind(p11.db,p11.dr)
    p10.d<- cbind(as.vector(u1*(1-u1)*SpDeriv1)*X,matrix(0,N,1))-p11.d
    p01.d<- cbind(as.vector(u2*(1-u2)*SpDeriv2)*X,matrix(0,N,1))-p11.d
    p00.d<- -p11.d-p10.d-p01.d
  }
  if(dim == p){ # dim == p
    p11.d<- p11.db
    p10.d<- as.vector(u1*(1-u1)*SpDeriv1)*X-p11.d
    p01.d<- as.vector(u2*(1-u2)*SpDeriv2)*X-p11.d
    p00.d<- -p11.d-p10.d-p01.d
  }
  return(list(p00.d=p00.d,p10.d=p10.d,p01.d=p01.d,p11.d=p11.d))
}