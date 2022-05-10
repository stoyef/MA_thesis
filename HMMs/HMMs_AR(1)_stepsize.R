# HMM WASP (HMM with autorrelation in state dependent process)

###
### AR(1) for step size
###

###
### read in data
#setwd(getSrcDirectory()[1]) # set wd to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # RStudio way

source("../autocorrelation/tsa_functions.R")
schwalben = read_schwalbe('../data/seeschwalbe/slick2-h1-h2.csv')

count=1
for (i in unique(schwalben$ID)){
  assign(paste("schwalbe_", count, sep=""), select_schwalbe(i, count))
  count=count+1
}

## work with schwalbe_77
library(moveHMM)
schwalbe_77 = prepData(schwalbe_77)


############## first: normal HMM using moveHMM, 2 states and 3 states
oldschool_2states <- fitHMM(schwalbe_77,2, c(20,30,2,4),angleDist = 'none')
oldschool_2states
plot(oldschool_2states)

oldschool_3states <- fitHMM(schwalbe_77,3, c(15,20,30,2,4,6),angleDist = 'none')
oldschool_3states
plot(oldschool_3states)

############# second: Include AR(1) in state dependent process
mllk_ar1_stepsize<-function(theta.star,x,N){
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  
  autocor <- plogis(theta.star[(N-1)*N+1:N])
  mu.step <- exp(theta.star[(N-1)*N+(N+1):(2*N)])
  sigma.step <- exp(theta.star[(N-1)*N+2*N+1:N])
  
  allprobs <- matrix(1,dim(x)[1],N)
  ind.step <- which(!is.na(x$step))[-c(1:1)] # change: we omit first step 
  # in order to always have the step in t-1

  for (j in 1:N){
    step.prob <- rep(1,dim(x)[1])
    # here comes the autocorrelation!
    mu.step_auto <- c(rep(1,1), # AR(1), more 1's before if higher order
                      (1-autocor[j])*mu.step[j] + autocor[j]*x$step[ind.step-1])
    step.prob[ind.step] <- dgamma(x$step[ind.step],
                                  shape=mu.step_auto[ind.step]^2/sigma.step[j]^2,
                                  scale=sigma.step[j]^2/mu.step_auto[ind.step]) # rephrase parameters
                      # here we have to choose mu.step_auto[ind.step], because
                      # we have an individual mu for each data point
    allprobs[,j] <- step.prob
  }
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:dim(x)[1]){
    foo <- phi%*%Gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}


###
### 2-state HMM
###
# starting values
theta = c(
  -2,-2, # values to construct TPM
  0.1,0.1, # autocorrelation [0,1]
  20,40, # means of step for each state [0,Inf)
  5, 5 # sd of step for each state [0,Inf)
) 
# transformation for unconstrained optimization
theta.star = c(
  theta[1:2], # values to construct TPM
  qlogis(theta[3:4]), # autocorrelation
  log(theta[5:6]), # step mean
  log(theta[7:8]) # step sd
)

mllk_ar1_stepsize(theta.star=theta.star, x=schwalbe_77[1:dim(schwalbe_77)[1]-1,], 2) # works

# minimize -logL
mod <- nlm(mllk_ar1_stepsize,theta.star,x=schwalbe_77[1:dim(schwalbe_77)[1]-1,],N=2,print.level=2,
           iterlim = 1000)
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



###
###
### state decoding
###
###
### global decoding, using viterbi from lecture
###

# We need TPM Gamma
autocor.viterbi <-function(x, Gamma, delta, autocor, 
                           mu.step, sigma.step){
  
  n <- dim(x)[1]
  allprobs <- matrix(1,n,2)
  ind <- which(!is.na(x$step))[-1] # change: we omit first step 
  # in order to always have the step in t-1
  
  allprobs[ind,] <- cbind(
    dgamma(x$step[ind],
           shape=mu.step[1]^2/sigma.step[1]^2,
           scale=sigma.step[1]^2/mu.step[1]),
    dgamma(x$step[ind],
           shape=mu.step[2]^2/sigma.step[2]^2,
           scale=sigma.step[2]^2/mu.step[2])
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
                                 Gamma, delta, autocor, mu.step, sigma.step)
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

hist(schwalbe_77$step, prob=T, breaks=40, xlab="Step size",
     ylim=c(0,0.8))

x <- seq(0,45,by=0.0005)
curve(delta[1]*dgamma(x, shape=mu.step[1]^2/sigma.step[1]^2,
                      scale=sigma.step[1]^2/mu.step[1]), 0,45, add=T,
      col=4,lwd=2)
curve(delta[2]*dgamma(x, shape=mu.step[2]^2/sigma.step[2]^2,
                      scale=sigma.step[2]^2/mu.step[2]), 0,45,add=T,
      col=7,lwd=2)
curve(delta[1]*dgamma(x, shape=mu.step[1]^2/sigma.step[1]^2,
                      scale=sigma.step[1]^2/mu.step[1])+
        delta[2]*dgamma(x, shape=mu.step[2]^2/sigma.step[2]^2,
                      scale=sigma.step[2]^2/mu.step[2]), 0,45,add=T,
      col=2,lwd=2)

# das sieht ja noch nicht so sinnvoll aus...





###
### 3-state HMM
###
# starting values
theta = c(
  rep(-2,6), # values to construct TPM
  0.1,0.2,0.2, # autocorrelation [0,1]
  20,25,35, # means of step for each state [0,Inf)
  6, 3, 2 # sd of step for each state [0,Inf)
) 
# transformation for unconstrained optimization
theta.star = c(
  theta[1:6], # values to construct TPM
  qlogis(theta[7:9]), # autocorrelation
  log(theta[10:12]), # step mean
  log(theta[13:15]) # step sd
)

mllk_ar1_stepsize(theta.star=theta.star, x=schwalbe_77, 3) # works

# minimize -logL
mod <- nlm(mllk_ar1_stepsize,theta.star,x=schwalbe_77[1:dim(schwalbe_77)[1]-1,],N=3,print.level=2,
           iterlim = 1000)
mod

## re-transformation to natural parameters
## TPM
N=3
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



###
###
### state decoding
###
###
### global decoding, using viterbi from lecture
###

# We need TPM Gamma
autocor.viterbi <-function(x, Gamma, delta, autocor, 
                           mu.step, sigma.step){
  
  n <- dim(x)[1]
  allprobs <- matrix(1,n,3)
  ind <- which(!is.na(x$step))[-c(1:1)] # change: we omit first step 
  # in order to always have the step in t-1
  
  allprobs[ind,] <- cbind(
    dgamma(x$step[ind],
           shape=mu.step[1]^2/sigma.step[1]^2,
           scale=sigma.step[1]^2/mu.step[1]),
    dgamma(x$step[ind],
           shape=mu.step[2]^2/sigma.step[2]^2,
           scale=sigma.step[2]^2/mu.step[2]),
    dgamma(x$step[ind],
           shape=mu.step[3]^2/sigma.step[3]^2,
           scale=sigma.step[3]^2/mu.step[3])
  )
  
  xi <- matrix(0,n,3)
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
                                 Gamma, delta, autocor, mu.step, sigma.step)
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

hist(schwalbe_77$step, prob=T, breaks=40, xlab="Step size",
     ylim = c(0,0.8))

x <- seq(0,45,by=0.0005)
curve(delta[1]*dgamma(x, shape=mu.step[1]^2/sigma.step[1]^2,
                      scale=sigma.step[1]^2/mu.step[1]), 0,45, add=T,
      col=4,lwd=2)
curve(delta[2]*dgamma(x, shape=mu.step[2]^2/sigma.step[2]^2,
                      scale=sigma.step[2]^2/mu.step[2]), 0,45,add=T,
      col=6,lwd=2)
curve(delta[3]*dgamma(x, shape=mu.step[3]^2/sigma.step[3]^2,
                      scale=sigma.step[3]^2/mu.step[3]), 0,45,add=T,
      col=7,lwd=2)
curve(delta[1]*dgamma(x, shape=mu.step[1]^2/sigma.step[1]^2,
                      scale=sigma.step[1]^2/mu.step[1])+
        delta[2]*dgamma(x, shape=mu.step[2]^2/sigma.step[2]^2,
                        scale=sigma.step[2]^2/mu.step[2])+
        delta[3]*dgamma(x, shape=mu.step[3]^2/sigma.step[3]^2,
                        scale=sigma.step[3]^2/mu.step[3]), 0,45, add=T,
      col=2,lwd=2)

# das sieht ja noch nicht so sinnvoll aus...





###
### next steps: gucken, warum, die einzelnen States so komisch sind...
###

