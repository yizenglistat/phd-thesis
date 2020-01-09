set.seed(1)
library(lsei)
library(splines2)

par(mfrow=c(1,1))
x = seq(-3,3,length.out=100)+rnorm(100,0,0.4)
y = seq(-3,3,length.out=100)^3
plot(x,y)
h = iSpline(x,knots=seq(-2.8,2.8,1),order=3,intercept=T)
#g = cbind(1,h[,-1])
g = cbind(1,h)
u = pnnls(g,y,1)$x
v = g%*%u
lines(x[order(x)],v[order(x)],lwd=2,col='red')
m = mSpline(x,knots=seq(-2.8,2.8,1),order=3,intercept = T)
#m = iSpline(x,knots=seq(-2.8,2.8,1),order=3,derivs = 1L)
n = m %*% u[-1]
lines(x[order(x)],n[order(x)],lwd=4,col='yellow')
lines(x[order(x)],3*(x[order(x)])^2,col='green')
