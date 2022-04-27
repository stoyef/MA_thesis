# HMM WASP (HMM with autorrelation in state dependent process)

###
### read in data
#setwd(getSrcDirectory()[1]) # set wd to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # RStudio way

source("../autocorrelation/tsa_functions.R")

schwalben = read_schwalbe('../data/seeschwalbe/slick2-h1-h2.csv')

head(schwalben)

count=1
for (i in unique(schwalben$ID)){
  assign(paste("schwalbe_", count, sep=""), select_schwalbe(i, count))
  count=count+1
}

## work with schwalbe_77
library(moveHMM)
schwalbe_77 = prepData(schwalbe_77)
head(schwalbe_77)


# compute -logL (from lecture code)

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
    allprobs[ind,j] <- dnorm(x[ind],mu[j],sigma[j]) # in this case Gaussian
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


###
### 1st aim: Clone model from Lawler et. al (2019)
###       (with exception of von Mises distribution instead of Wrappen Cauchy)
###
### (Code adapted from HMM lecture)

mllk_ar1<-function(theta.star,x,N){
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  
  autocor <- plogis(theta.star[(N-1)*N+1:N])
  mu.step <- exp(theta.star[(N-1)*N+(N+1):(2*N)])
  sigma.step <- exp(theta.star[(N-1)*N+2*N+1:N])
  
  # same transformation as in function w2n of moveHMM package
  mu.angle <- Arg(theta.star[(N-1)*N+3*N+1:N]+1i*theta.star[(N-1)*N+4*N+1:N])
  kappa.angle <- sqrt(theta.star[(N-1)*N+3*N+1:N]^2+theta.star[(N-1)*N+4*N+1:N]^2)
  
  allprobs <- matrix(1,dim(x)[1],N)
  ind.step <- which(!is.na(x$step))[-1] # change: we omit first step 
  # in order to always have the step in t-1
  ind.angle <- which(!is.na(x$angle))
  
  for (j in 1:N){
    step.prob <- rep(1,dim(x)[1])
    angle.prob <- rep(1,dim(x)[1]) # missing observations stay at value 1
    # here comes the autocorrelation!
    mu.step_auto <- (1-autocor[j])*mu.step[j] + autocor[j]*x$step[ind.step-1]
    step.prob[ind.step] <- dgamma(x$step[ind.step],
                        shape=mu.step_auto[j]^2/sigma.step[j]^2,
                        scale=sigma.step[j]^2/mu.step_auto[j]) # rephrase parameters
    angle.prob[ind.angle] <- dvm(x$angle[ind.angle],mu.angle[j],kappa.angle[j]) 
    allprobs[,j] <- step.prob * angle.prob
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


# starting values
theta = c(
  -2,-2, # values to construct TPM
  0.3,0.5, # autocorrelation [0,1]
  20,35, # means of step for each state [0,Inf)
  5, 10, # sd of step for each state [0,Inf)
  0, 0, # means of angle for each state [-pi,pi]
  2, 4 # kappa of angle for each state [0,Inf)
) 
# transformation for unconstrained optimization
theta.star = c(
  theta[1:2], # values to construct TPM
  qlogis(theta[3:4]), # autocorrelation
  log(theta[5:6]), # step mean
  theta[7:8], # step sd
  # for parameter of von Mises distribution: same transformation as in 
  # function n2w of moveHMM package
  theta[9:10] * cos(theta[11:12]), # angle mean
  theta[9:10] * sin(theta[11:12]) # angle kappa
)

mllk_ar1(theta.star=theta.star, x=schwalbe_77, 2) # works

# minimize -logL
mod <- nlm(mllk_ar1,theta.star,x=schwalbe_77[2:4480,],N=2,print.level=2)
mod

## re-transformation to natural parameters
## TPM
N=2
Gamma <- diag(N)
Gamma[!Gamma] <- exp(mod$estimate[1:((N-1)*N)])
Gamma <- Gamma/rowSums(Gamma)
Gamma

delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
delta

# autocorrelation
plogis(mod$estimate[(N-1)*N+1:N])

# step length
mu.step <- exp(mod$estimate[(N-1)*N+(N+1):(2*N)])
mu.step
sigma.step <- exp(mod$estimate[(N-1)*N+2*N+1:N])
sigma.step

# turning angles
mu.angle <- Arg(mod$estimate[(N-1)*N+3*N+1:N]+1i*mod$estimate[(N-1)*N+4*N+1:N])
mu.angle
kappa.angle <- sqrt(mod$estimate[(N-1)*N+3*N+1:N]^2+mod$estimate[(N-1)*N+4*N+1:N]^2)
kappa.angle





###
###
### ok, this seems to kinda work in principle but seems very unstable
### also, the estimated autocorrelation behaves very differently, depending
### on the starting values
###
###

### next try: fix mean turning angle to 0 to reduce parameters


mllk_ar1<-function(theta.star,x,N){
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  
  autocor <- plogis(theta.star[(N-1)*N+1:N])
  mu.step <- exp(theta.star[(N-1)*N+(N+1):(2*N)])
  sigma.step <- exp(theta.star[(N-1)*N+2*N+1:N])
  
  # change to before: mu.angle fixed
  mu.angle <- rep(0,N)
  kappa.angle <- exp(theta.star[(N-1)*N+3*N+1:N])
  
  allprobs <- matrix(1,dim(x)[1],N)
  ind.step <- which(!is.na(x$step))[-1] # change: we omit first step 
  # in order to always have the step in t-1
  ind.angle <- which(!is.na(x$angle))
  
  for (j in 1:N){
    step.prob <- rep(1,dim(x)[1])
    angle.prob <- rep(1,dim(x)[1]) # missing observations stay at value 1
    # here comes the autocorrelation!
    mu.step_auto <- (1-autocor[j])*mu.step[j] + autocor[j]*x$step[ind.step-1]
    step.prob[ind.step] <- dgamma(x$step[ind.step],
                                  shape=mu.step_auto[j]^2/sigma.step[j]^2,
                                  scale=sigma.step[j]^2/mu.step_auto[j]) # rephrase parameters
    angle.prob[ind.angle] <- dvm(x$angle[ind.angle],mu.angle[j],kappa.angle[j]) 
    allprobs[,j] <- step.prob * angle.prob
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


# starting values
theta = c(
  -2,-2, # values to construct TPM
  0.3,0.5, # autocorrelation [0,1]
  20,35, # means of step for each state [0,Inf)
  5, 10, # sd of step for each state [0,Inf)
  #0, 0, # means of angle for each state [-pi,pi]
  2, 4 # kappa of angle for each state [0,Inf)
) 
# transformation for unconstrained optimization
theta.star = c(
  theta[1:2], # values to construct TPM
  qlogis(theta[3:4]), # autocorrelation
  log(theta[5:6]), # step mean
  theta[7:8], # step sd
  # for parameter of von Mises distribution: same transformation as in 
  # function n2w of moveHMM package
  #theta[9:10] * cos(theta[11:12]), # angle mean
  log(theta[9:10]) # angle kappa
)

mllk_ar1(theta.star=theta.star, x=schwalbe_77, 2) # works

# minimize -logL
mod <- nlm(mllk_ar1,theta.star,x=schwalbe_77[2:4480,],N=2,print.level=2)
mod

## re-transformation to natural parameters
## TPM
N=2
Gamma <- diag(N)
Gamma[!Gamma] <- exp(mod$estimate[1:((N-1)*N)])
Gamma <- Gamma/rowSums(Gamma)
Gamma

delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
delta

# autocorrelation
autocor <- plogis(mod$estimate[(N-1)*N+1:N])
autocor

# step length
mu.step <- exp(mod$estimate[(N-1)*N+(N+1):(2*N)])
mu.step
sigma.step <- exp(mod$estimate[(N-1)*N+2*N+1:N])
sigma.step

# turning angles
kappa.angle <- exp(mod$estimate[(N-1)*N+3*N+1:N])
kappa.angle



###
###
### state decoding
###
###
### global decoding, using viterbi from lecture
###

# We need TPM Gamma
autocor.viterbi <-function(x, Gamma, delta, autocor, 
                           mu.step, sigma.step, kappa.angle){
 
  n <- dim(x)[1]
  allprobs <- matrix(1,n,2)
  ind <- which(!is.na(x$step))[-1] # change: we omit first step 
  # in order to always have the step in t-1
  
  allprobs[ind,] <- cbind(
    dgamma(x$step[ind],
           shape=mu.step[1]^2/sigma.step[1]^2,
           scale=sigma.step[1]^2/mu.step[1]) *
      dvm(x$angle[ind],0,kappa.angle[1]), # mu.angle fixed as 0
    dgamma(x$step[ind],
           shape=mu.step[2]^2/sigma.step[2]^2,
           scale=sigma.step[2]^2/mu.step[2]) *
      dvm(x$angle[ind],0,kappa.angle[2]) # mu.angle fixed as 0
  )
  
  xi <- matrix(0,n,2)
  foo <- delta*allprobs[1,]
  xi[1,] <- foo/sum(foo)
  
  for (t in 2:n){
    foo <- apply(xi[t-1,]*Gamma, 2, max) * allprobs[t,]
    xi[t,] <- foo/sum(foo)
  }
  
  iv <- numeric(n)
  iv[n] <- which.max(xi[n,])
  
  for (t in (n-1):1){
    iv[t] <- which.max(Gamma[,iv[t+1]] * xi[t,])
  }
  return(iv)
}

states_global <- autocor.viterbi(schwalbe_77, 
                                 Gamma, delta, autocor, mu.step, sigma.step,
                                 kappa.angle)
states_global


###
###
###
### Visualization
###
###
###


# Visualization of the different distributions

# step
plot(schwalbe_77$step, type='l')
points(schwalbe_77$step, pch=19, col=states_global+1) # vernünftige viz fehlt

hist(schwalbe_77$step, prob=T, breaks=20, xlab="Step size")

x <- seq(0,45,by=0.0005)
curve(delta[1]*dgamma(x, shape=mu.step[1]^2/sigma.step[1]^2,
                    scale=sigma.step[1]^2/mu.step[1]), 0,45, add=T)
curve(delta[2]*dgamma(x, shape=mu.step[2]^2/sigma.step[2]^2,
             scale=sigma.step[2]^2/mu.step[2]), 0,45,add=T)

# das sieht ja noch nicht so sinnvoll aus...


# angle
plot(schwalbe_77$angle, type='l')
points(schwalbe_77$angle, pch=19, col=states_global+1) # vernünftige viz fehlt

hist(schwalbe_77$angle, prob=T, breaks=70, xlab="angle",xlim=c(-1,1))

x <- seq(-1,1,by=0.005)
curve(delta[1]*dvm(x, 0,kappa.angle[1]), -1,1,add=T)
curve(delta[2]*dvm(x ,0,kappa.angle[2]), -1,1,add=T)
curve(delta[2]*dvm(x ,0,kappa.angle[2]), -1,1)

## ach manno, irgendwas ist da Quark, nochmal alles durchgehen als nächstes :)





