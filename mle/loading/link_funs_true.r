# links is a list collection of user-defined link functions

# monotone
eta1=function(u){
  return(1.4*u-9)
}

eta2=function(u){
 #return(1-3/pcauchy(u))
 return(-10+u^3/25+u^2/18+u/18)
}

g1<-function(u){return(1/(1+exp(-eta1(u))))}
g2<-function(u){return(1/(1+exp(-eta2(u))))}

#g1<-function(u){return(1/(1+exp(4-u/(1+abs(u)))))}

#g2<-function(u){return(1/(1+exp(4-atan(u))))}

#g2<-function(u){return(1/(1+exp(3-u^3/25-u^2/18-u/18)))}
#g2<- function(u){return(0.95/(1+exp(-5*u+20)) + 0.05/(1+exp(-5*u-10)))}
#g2<- function(u){return(0.25/(1+exp(-50*u-15)) + 0.75/(1+exp(-10*u+15)))}
#g3<-function(u){return((sin(pi*(u-0.3))+1.3)/(20+30*(u-0.3)^2*(sign(u-0.3)+1)))}

links <- list(g1,g2)