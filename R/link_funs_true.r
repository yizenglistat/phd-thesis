# links is a list collection of user-defined link functions

# monotone
# eta1=function(u){
#   return(1*sign(u)*(abs(u))^(1.2)-8.5)
# }

# eta2=function(u){
#  #return()
#  return(0.5*sign(u)*(abs(u))^(1.4)-6)
# }

#g1<-function(u){return(1/(1+exp(-eta1(u))))}
#g2<-function(u){return(1/(1+exp(-eta2(u))))}



g1 <- function(u) 1/(1+exp(55-40*atan(u)))#pnorm(u,8,3)
g2 <- function(u) pcauchy(u,7,2/3)

eta1 <- function(u) log(g1(u)/(1-g1(u)))
eta2 <- function(u) log(g2(u)/(1-g2(u)))

etas <- list(eta1,eta2)

#g1<-function(u){return(1/(1+exp(4-u/(1+abs(u)))))}

#g2<-function(u){return(1/(1+exp(4-atan(u))))}

#g2<-function(u){return(1/(1+exp(3-u^3/25-u^2/18-u/18)))}
#g2<- function(u){return(0.95/(1+exp(-5*u+20)) + 0.05/(1+exp(-5*u-10)))}
#g2<- function(u){return(0.25/(1+exp(-50*u-15)) + 0.75/(1+exp(-10*u+15)))}
#g3<-function(u){return((sin(pi*(u-0.3))+1.3)/(20+30*(u-0.3)^2*(sign(u-0.3)+1)))}

links <- list(g1,g2)