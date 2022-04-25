
#### This is code from the lecture Hidden Markov Models ####
#### by Roland Langrock (March 2021) #######################




### Old Faithful mixture model 
OF<-read.table("http://www.rolandlangrock.com//Misc//OF.dat") 

# Likelihood of mixture model
l<-function(theta,x){ 
  mu<-theta[1:2] 
  sigma<-theta[3:4] 
  pi<-theta[5] 
  logl<-sum(log(pi*dnorm(x,mu[1],sigma[1])+ 
                  (1-pi)*dnorm(x,mu[2],sigma[2]))) 
  return(-logl) 
} 
mod<-nlminb(c(2,4.5,1,2,0.3),l,x=OF$eruptions, 
            lower=c(0,0,0,0,0),upper=c(Inf,Inf,Inf,Inf,1)) 
mod



### delta 
Gamma<-matrix(c(0.65,0.35,0.21,0.79),byrow=TRUE,nrow=2) 
delta <- solve(t(diag(2)-Gamma+1),c(1,1)) 
delta 

delta<-c(1,0) 
for (k in 1:30){ 
  delta<-delta%*%Gamma 
  print(delta) 
} 

g1<-matrix(c(0.1,0.2,0.3,0.4,0.1,0.6,0.2,0.1,0.2,0.2,0.4,0.2,0.3,0.2,0.2,0.3),
           byrow = TRUE, nrow=4)
d1<-solve(t(diag(4)-g1+1), c(1,1,1,1))
d1
d1%*%g1


### simulate from HMM 
n<-50 
x<-s<-rep(NA,n) 
Gamma<-matrix(c(0.9,0.1,0.1,0.9),nrow=2) 
delta<-c(0.5,0.5) 
mu<-c(5,14) 
sigma<-c(2,3) 
s[1]<-sample(1:2,size=1,prob=delta) 
x[1]<-rnorm(1,mu[s[1]],sigma[s[1]]) 
for (t in 2:50){ 
  s[t]<-sample(1:2,size=1,prob=Gamma[s[t-1],]) 
  x[t]<-rnorm(1,mu[s[t]],sigma[s[t]]) 
}


### likelihood evaluation (2-state Poisson HMM)

L<-function(theta,x){
  Gamma <- diag(theta[1:2])
  Gamma[1,2] <- 1-Gamma[1,1]
  Gamma[2,1] <- 1-Gamma[2,2] 
  delta <- solve(t(diag(2)-Gamma+1),c(1,1))
  lambda <- theta[3:4]
  allprobs <- cbind(dpois(x,lambda[1]),dpois(x,lambda[2]))
  foo <- delta%*%diag(allprobs[1,])
  for (t in 2:length(x)){
    foo <- foo%*%Gamma%*%diag(allprobs[t,])
  }
  return(sum(foo))
}

quakes<-read.table("http://www.rolandlangrock.com/Misc/earthquakes.txt",header=TRUE)
theta<-c(0.85,0.81,15,25) # a fairly randomly chosen parameter vector
L(theta,x=quakes$count) # likelihood
log(L(theta,x=quakes$count)) # log-likelihood



### fitting a 2-state Poisson HMM to the earthquake data

mllk<-function(theta.star,x){
  theta <- c(plogis(theta.star[1]),plogis(theta.star[2]),
             exp(theta.star[3]),exp(theta.star[4]))
  Gamma <- diag(theta[1:2])
  Gamma[1,2] <- 1-Gamma[1,1]
  Gamma[2,1] <- 1-Gamma[2,2] 
  delta <- solve(t(diag(2)-Gamma+1),c(1,1))
  lambda <- theta[3:4]
  allprobs <- cbind(dpois(x,lambda[1]),dpois(x,lambda[2]))
  foo <- delta%*%diag(allprobs[1,])
  for (t in 2:length(x)){
    foo <- foo%*%Gamma%*%diag(allprobs[t,])
  }
  return(-log(sum(foo))) # returns -log(L) since nlm() can only minimise!
}


quakes<-read.table("http://www.rolandlangrock.com/Misc/earthquakes.txt",header=TRUE) 
dgammas<-c(0.8,0.8) # starting values for diagonal entries in t.p.m.
lambdas<-c(15,25) # starting values for the state-dependent lambdas
theta.star <- c(qlogis(dgammas),log(lambdas))
s<-Sys.time()
mod<-nlm(mllk,theta.star,x=quakes$count,print.level=2)
Sys.time()-s

c(plogis(mod$estimate[1:2]),exp(mod$estimate[3:4]))



### same but now addressing potential underflow

mllk<-function(theta.star,x){
  theta <- c(plogis(theta.star[1]),plogis(theta.star[2]),
             exp(theta.star[3]),exp(theta.star[4]))
  Gamma <- diag(theta[1:2])
  Gamma[1,2] <- 1-Gamma[1,1]
  Gamma[2,1] <- 1-Gamma[2,2] 
  delta <- solve(t(diag(2)-Gamma+1),c(1,1))
  lambda <- theta[3:4]
  allprobs <- cbind(dpois(x,lambda[1]),dpois(x,lambda[2]))
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:length(x)){
    foo <- phi%*%Gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}

theta.star<-c(qlogis(c(0.8,0.8)),log(c(15,25))) 
nlm(mllk,theta.star,x=quakes$count,print.level=2)



### local maxima check 

llks <- rep(NA,100)
mods <- vector("list")
for (k in 1:100){
  theta.star <- c(qlogis(runif(2,0,1)),log(runif(2,10,30)))
  mods[[k]] <- nlm(mllk,theta.star,x=quakes$count,stepmax=5)
  llks[k] <- -mods[[k]]$minimum
}

llks



### fitting a 2-state gamma HMM to the muskox data

mllk<-function(theta.star,x){
  theta <- c(plogis(theta.star[1]),plogis(theta.star[2]),
             exp(theta.star[3]),exp(theta.star[4]),
             exp(theta.star[5]),exp(theta.star[6]))
  Gamma <- diag(theta[1:2])
  Gamma[1,2] <- 1-Gamma[1,1]
  Gamma[2,1] <- 1-Gamma[2,2] 
  delta <- solve(t(diag(2)-Gamma+1),c(1,1))
  mu <- theta[3:4]
  sigma <- theta[5:6]
  allprobs <- matrix(1,length(x),2)
  ind<-which(!is.na(x))
  allprobs[ind,] <- cbind(
    dgamma(x[ind],shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1]),
    dgamma(x[ind],shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2]))
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:length(x)){
    foo <- phi%*%Gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}

muskox <-read.csv("http://www.rolandlangrock.com//Misc//muskox.csv")
library(moveHMM) # the package moveHMM needs to be installed first!
data<-prepData(muskox,type="UTM")
mod



### local maxima check

llks<-rep(NA,100)
mods<-vector("list")
for (k in 1:100){
  theta.star <- c(qlogis(runif(2,0,1)),log(runif(2,0,1000)),log(runif(2,0,1000))) 
  mods[[k]] <- nlm(mllk,theta.star,x=data$step,stepmax=5)
  llks[k] <- -mods[[k]]$minimum
}



### confidence intervals
theta.star <- c(qlogis(runif(2,0,1)),log(runif(2,0,1000)),log(runif(2,0,1000))) 
mod<-nlm(mllk,theta.star,x=data$step,hessian=TRUE)
sds<-sqrt(diag(solve(mod$hessian)))
CIs<-cbind(mod$estimate-1.96*sds,mod$estimate+1.96*sds)
CIs

rbind(plogis(CIs[1,]),plogis(CIs[2,]),exp(CIs[3,]),exp(CIs[4,]),exp(CIs[5,]),exp(CIs[6,]))



### fitting N-state gamma HMMs

mllk<-function(theta.star,x,N){
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  mu <- exp(theta.star[(N-1)*N+1:N])
  sigma <- exp(theta.star[(N-1)*N+(N+1):(2*N)])
  allprobs <- matrix(1,length(x),N)
  ind <- which(!is.na(x))
  for (j in 1:N){
    allprobs[ind,j] <- dgamma(x[ind],
                              shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j])
  }
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:length(x)){
    foo <- phi%*%Gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}

N=3
theta.star<-c(rep(-2,(N-1)*N),log(seq(5,500,length=N)),log(seq(5,500,length=N))) 
mod<-nlm(mllk,theta.star,x=data$step,N=N,print.level=2,iterlim=10000)




### live coding with caracaras data

# caracaras (that's a bird!) acceleration data (the higher the values, the more active is the bird - highest values indicate flight mode)
cara <- read.csv("http://www.rolandlangrock.com//Misc//caracaras.csv")$x

mllk<-function(theta.star,x,N){
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  mu <- theta.star[(N-1)*N+1:N]
  sigma <- exp(theta.star[(N-1)*N+(N+1):(2*N)])
  allprobs <- matrix(1,length(x),N)
  ind <- which(!is.na(x))
  for (j in 1:N){
    allprobs[ind,j] <- dnorm(x[ind],mu[j],sigma[j])
  }
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:length(x)){
    foo <- phi%*%Gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}

N=3
theta.star<-c(rep(-2,(N-1)*N),-5,-3.8,-1.3,log(c(1,2,1))) 
mod<-nlm(mllk,theta.star,x=log(cara),N=N,print.level=2,iterlim=10000,stepmax=5)

# estimated unconstrained parameter vector:
theta.star <- mod$estimate
theta.star

# transform theta.star back to the constrained parameter space:
Gamma <- diag(N)
Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
Gamma <- Gamma/rowSums(Gamma)
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
mu <- theta.star[(N-1)*N+1:N]
sigma <- exp(theta.star[(N-1)*N+(N+1):(2*N)])
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))

# display estimated parameters:
round(Gamma,3)
round(mu,3)
round(sigma,3)

hist(log(cara),breaks=50,col="lightgrey",main="",xlab="log(ODBA)",ylab="density",prob=TRUE)
for (j in 1:3){
  curve(delta[j]*dnorm(x,mu[j],sigma[j]),from=-7,to=2,add=TRUE,lwd=2,col=j)
}




### fitting 2-8 state Gaussian HMMs to the caracaras data (log(ODBA) values)

cara <- read.csv("http://www.rolandlangrock.com//Misc//caracaras.csv")$x

mllk<-function(theta.star,x,N){
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  mu <- theta.star[(N-1)*N+1:N]
  sigma <- exp(theta.star[(N-1)*N+(N+1):(2*N)])
  allprobs <- matrix(1,length(x),N)
  ind <- which(!is.na(x))
  for (j in 1:N){
    allprobs[ind,j] <- dnorm(x[ind],mu[j],sigma[j])
  }
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:length(x)){
    foo <- phi%*%Gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}

allmods<-vector("list")
AIC<-vector("list")

for (N in 2:8){
  theta.star<-c(rep(-2,(N-1)*N),seq(-4.5,-1.6,length=N),log(rep(2,N))) 
  allmods[[N]]<-nlm(mllk,theta.star,x=log(cara),N=N,print.level=2,iterlim=10000)
  AIC[[N]]<-2*allmods[[N]]$minimum+2*length(allmods[[N]]$estimate)
}

cbPalette <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",1)
for (N in 8:2){
  theta.star <- allmods[[N]]$estimate
  theta.star
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  mu <- theta.star[(N-1)*N+1:N]
  sigma <- exp(theta.star[(N-1)*N+(N+1):(2*N)])
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  
  #windows()
  
  hist(log(cara),breaks=70,col="lightgrey",main=paste("N=",N),xlab="log(ODBA)",ylab="density",prob=TRUE)
  for (j in 1:N){
    curve(delta[j]*dnorm(x,mu[j],sigma[j]),from=-7,to=2,add=TRUE,lwd=2,col=cbPalette[j],n=1000)
  }
  
  text(-2,0.6,paste("AIC=",round(AIC[[N]],1)),cex=2)
}



### Viterbi for the 4-state Gaussian HMM fitted to the caracaras data

viterbi<-function(x,mu,sigma,Gamma,delta,N){
  n <- length(x)
  allprobs <- matrix(1,n,N)
  ind <- which(!is.na(x))
  for (j in 1:N){
    allprobs[ind,j] <- dnorm(x[ind],mu[j],sigma[j])
  } 
  xi <- matrix(0,n,N)
  foo <- delta*allprobs[1,]
  xi[1,] <- foo/sum(foo)
  for (t in 2:n){
    foo <- apply(xi[t-1,]*Gamma,2,max)*allprobs[t,]
    xi[t,] <- foo/sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n,])
  for (t in (n-1):1){
    iv[t] <- which.max(Gamma[,iv[t+1]]*xi[t,])
  }
  iv
}

N=3
theta.star <- allmods[[N]]$estimate
theta.star
Gamma <- diag(N)
Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
Gamma <- Gamma/rowSums(Gamma)
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
mu <- theta.star[(N-1)*N+1:N]
sigma <- exp(theta.star[(N-1)*N+(N+1):(2*N)])
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
states <- viterbi(log(cara),mu,sigma,Gamma,delta,N)

plot(log(cara),pch=19,col=cbPalette[states])


#### Chapter 7 #### -> state-space models

mllk <- function(theta.star,x,m,bm){
  phi <- plogis(theta.star[1])
  sigma <- exp(theta.star[2])
  beta <- exp(theta.star[3])
  b <- seq(-bm,bm,length=m+1)          # specify boundaries of m intervals
  h <- b[2]-b[1]                       # h is the length of each interval
  bstar <- (b[-1]+b[-(m+1)])*0.5       # midpoints of the m intervals
  Gamma <- matrix(0,m,m) 
  for (i in 1:m){
    Gamma[i,] <- h*dnorm(bstar,phi*bstar[i],sigma) # m*m t.p.m. of the approx. HMM
  }
  delta <- h*dnorm(bstar,0,sigma/sqrt(1-phi^2))   # stat. initial distribution
  foo <- delta*dpois(x[1],exp(bstar)*beta)
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:length(x)){
    foo <- phi%*%Gamma*dpois(x[t],exp(bstar)*beta)
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}

quakes<-read.table("http://www.rolandlangrock.com/Misc/earthquakes.txt",header=TRUE)

theta.star <- c(qlogis(0.8),log(0.2),log(20))
s=Sys.time()
mod<-nlm(mllk,theta.star,x=quakes$count,m=150,bm=1.5,print.level=2)
Sys.time()-s

c(plogis(mod$estimate[1]),exp(mod$estimate[2]),exp(mod$estimate[3]))

