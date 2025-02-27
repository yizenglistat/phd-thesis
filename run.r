# r script
rm(list=ls(all=TRUE)); source('./R/loading.r')
seed_number = 2000; set.seed(seed_number)

N        <- 3000   # sample size
c        <- 5      # group size
ord      <- 4      # order of splines
niknots  <- 10     # number of interior knots
Se       <- c(0.95,0.95) # sensitivity
Sp       <- c(0.95,0.95) # specificity
X  <- cbind(rnorm(N),
            rnorm(N),
            rbinom(N,1,0.1),
            rbinom(N,1,0.2),
            rbinom(N,1,0.3),
            rbinom(N,1,0.5)) # covariate
links<-list(links[[1]],links[[2]]) # true link functions g_k's

# true responese and data generated
Y <- simulate_y(X,beta,links,delta) # need true Y status to simulate data
cj<-sample(c(c(rep(4,125),rep(5,500))))

DATA <- simulate_data(X,Y,N,cj,Se,Sp) # summary of simulated data
cj <- DATA$cj # group setting
data <- DATA$data
colSums(Y)/N

out_list<-mle(X, cj, data, Se, Sp, ord, niknots, verbose=TRUE, seed=seed_number)
#*******************************************************************#
#                        Simulation Graphing
#*******************************************************************#
par(mfrow=c(2,2))
u <- sort(X%*%c(1,beta))
alpha1 <- pnnls(Bsp(u,ord,niknots),eta1(u),1)$x
alpha2 <- pnnls(Bsp(u,ord,niknots),eta2(u),1)$x
#alpha1 <- out$alpha1
#alpha2 <- out$alpha2

plot(u,g(u,alpha1,ord,niknots),'l',ylim=c(0,1),col='blue',lty=2)
lines(u,g1(u), type='l', col='green',lty=2)
plot(u,g(u,alpha2,ord,niknots),'l',ylim=c(0,1),col='blue',lty=2)
lines(u,g2(u), type='l', col='green',lty=2)


plot(u,eta.Bsp(u,alpha1,ord,niknots),'l',col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5),lty=4,lwd=2)
lines(u,eta1(u), type='l', col=rgb(red = 0, green = 1, blue = 0, alpha = 0.5),lty=4,lwd=2)
plot(u,eta.Bsp(u,alpha2,ord,niknots),'l',col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5),lty=4,lwd=2)
lines(u,eta2(u), type='l', col=rgb(red = 0, green = 1, blue = 0, alpha = 0.5),lty=4,lwd=2)

#fig_data = data.frame(u=c(u,u),g=c(g1(u),g2(u)),type=c(rep('g1',length(u)),rep('g2',length(u))))

#xyplot(g ~ u | type, fig_data, type = "l", pch=20, layout=c(2,1))

output_figure(c(alpha1,alpha2), beta, ord, niknots, etas, plot='eta')
#*******************************************************************#

#*******************************************************************#
#                          Data Application
#*******************************************************************#
# female_swab <- read.csv("./female_swab.csv")
# female_swab <- female_swab[female_swab$cj==4,]
# X <- as.matrix(female_swab[c("Age","Race","Risk.New.Partner",
#                              "Risk.Multiple.Partners",
#                              "Risk.Contact",
#                              "Symptoms")])
# cj<-c(as.numeric(table(female_swab$Pool.ID)),
#       rep(1,length(female_swab$Pool.ID[is.na(female_swab$Pool.ID)])))
# data <- female_swab[c("CT.Result","GC.Result","CT.Group","GC.Group")]
# N <- dim(X)[1]
# r <- 4 # order of splines
# m <- 4 # number of interior knots
# Se <- c(0.992,0.992)
# Sp <- c(0.992,0.992)
# Y <- SIM_true(X,beta,links,delta) # need true Y status to simulate data
# colSums(Y)/N

# female_urine <- read.csv("./female_urine.csv")
# X <- as.matrix(female_urine[c("Age","Race","Risk.New.Partner",
#                              "Risk.Multiple.Partners",
#                              "Risk.Contact",
#                              "Symptoms")])
# cj<-c(as.numeric(table(female_urine$Pool.ID)),
#       rep(1,length(female_urine$Pool.ID[is.na(female_urine$Pool.ID)])))
# data <- female_urine[c("CT.Result","GC.Result","CT.Group","GC.Group")]
# N <- dim(X)[1]
# r <- 4 # order of splines
# m <- 4 # number of interior knots
# Se <- c(0.947,0.913)
# Sp <- c(0.989,0.993)

#female_urine <- read.csv("./female_urine.csv")

####################### comment
X<-as.matrix(read.csv("./X.csv")[,-1])
data<-as.matrix(read.csv("./data.csv")[,-1])
cj<-read.csv("./cj.csv")[,-1]
age<-X[,1]

Agesq<-X[,1]^2
X<-cbind(Agesq,X)
X[,1]<-(X[,1]-mean(X[,1]))/(max(X[,1])-min(X[,1]))
X[,2]<-(X[,2]-mean(X[,2]))/(max(X[,2])-min(X[,2]))
X<-X[,-9]
X[,3:10]<-ifelse(X[,3:10]==1,0.5,-0.5)
X <- X[,c(2,1,3:10)]

N <- dim(X)[1]
r <- 3 # order of splines
m <- 4 # number of interior knots
Se <- c(0.942,0.992)
Sp <- c(0.976,0.987)


#*******************************************************************#
#                             Graphing
#*******************************************************************#
final.par<-mle(X, cj, data, Se, Sp, r, m, eps=5e-3, verbose=T, maxiter=20)
out<-capture.output(final.par<-mle(X, cj, data, Se, Sp, r, m, eps=5e-3, verbose=T, maxiter=20))
par<-capture.output(final.par)
cat("basis splines with age square", out, file="out_ms.txt", sep="\n", append=T)
# 
# 
# out<-capture.output(final.par<-MLE(X, cj, data, Se, Sp, r, m, eps=5e-3, verbose=T, maxiter=15))
# par<-capture.output(final.par)
# cat("monotone splines without age square", out, file="out_ms.txt", sep="\n", append=TRUE)
# cat("monotone splines without age square", par, file="out_ms.txt", sep="\n", append=TRUE)
# 


# graph the estimated link functions
par(mfrow=c(1,2))
est.beta<-final.par$final.beta

#u<-X%*%c(1,est.beta)
sort.u<-sort(X%*%est.beta)
#u<-u/sqrt(sum(u^2))
#sort.u<-sort(X%*%c(1,est.beta))

sol1<-final.par$final.alpha1
sol2<-final.par$final.alpha2
 
plot(sort.u,g(sort.u,sol1,r,m),type='l',col='red',lwd=2,
     xlab=expression(u),
     ylab=expression(hat(g)[1](u)));
plot(sort.u,g(sort.u,sol2,r,m),type='l',col='red',lwd=2,
     xlab=expression(u),
     ylab=expression(hat(g)[2](u)));

par(mfrow=c(1,2))
plot(sort.u,log(g(sort.u,sol1,r,m)/(1-g(sort.u,sol1,r,m))),type='l',col='red',lwd=2,
     xlab=expression(u),
     ylab=expression(hat(g)[1](u)));
plot(sort.u,log(g(sort.u,sol2,r,m)/(1-g(sort.u,sol2,r,m))),type='l',col='red',lwd=2,
     xlab=expression(u),
     ylab=expression(hat(g)[2](u)));

