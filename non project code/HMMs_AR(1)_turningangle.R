# HMM WASP (HMM with autorrelation in state dependent process)

###
### AR(1) for turning angle
###


###
### read in data
#setwd(getSrcDirectory()[1]) # set wd to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # RStudio way

source("../autocorrelation/tsa_functions.R")
schwalben = read_schwalbe('../datasets/seeschwalbe/slick2-h1-h2.csv')

count=1
for (i in unique(schwalben$ID)){
  assign(paste("schwalbe_", count, sep=""), select_schwalbe(i, count))
  count=count+1
}

## work with schwalbe_77
library(moveHMM)
schwalbe_77 = prepData(schwalbe_77)



############# Include AR(1) in state dependent process
mllk_ar1_turningangle<-function(theta.star,x,N){
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  
  autocor <- plogis(theta.star[(N-1)*N+1:N])
  
  # same transformation as in function w2n of moveHMM package
  mu.angle <- Arg(theta.star[(N-1)*N+(N+1):(2*N)]+1i*theta.star[(N-1)*N+2*N+1:N])
  kappa.angle <- sqrt(theta.star[(N-1)*N+(N+1):(2*N)]^2+theta.star[(N-1)*N+2*N+1:N]^2)
  
  allprobs <- matrix(1,dim(x)[1],N)
  ind.angle <- which(!is.na(x$angle))[-c(1:1)] # change: we omit first step 
  # in order to always have the step in t-1
  
  for (j in 1:N){
    angle.prob <- rep(1,dim(x)[1])
    # here comes the autocorrelation!
    mu.angle_auto <- c(rep(NA,1), # AR(1), more 1's before if higher order
                      (1-autocor[j])*mu.angle[j] + autocor[j]*x$angle[ind.angle-1])
    angle.prob[ind.angle] <- dvm(x$angle[ind.angle],mu.angle_auto[ind.angle],kappa.angle[j])
                                  
    # here we have to choose mu.step_auto[ind.step], because
    # we have an individual mu for each data point
    allprobs[,j] <- angle.prob
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
  0.1,-0.1, # means of angle for each state [-pi,pi]
  300, 1000 # kappa of step for each state [0,Inf)
) 
# transformation for unconstrained optimization
theta.star = c(
  theta[1:2], # values to construct TPM
  qlogis(theta[3:4]), # autocorrelation
  theta[5:6] * cos(theta[7:8]), # angle mean
  theta[5:6] * sin(theta[7:8]) # angle kappa
)

mllk_ar1_turningangle(theta.star=theta.star, x=schwalbe_77[2:dim(schwalbe_77)[1],], 2) # works

# minimize -logL
mod <- nlm(mllk_ar1_turningangle,theta.star,x=schwalbe_77[2:dim(schwalbe_77)[1],],N=2,print.level=2,
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

# angle
mu.angle <- Arg(mod$estimate[(N-1)*N+(N+1):(2*N)]+1i*mod$estimate[(N-1)*N+2*N+1:N])
kappa.angle <- sqrt(mod$estimate[(N-1)*N+(N+1):(2*N)]^2+mod$estimate[(N-1)*N+2*N+1:N]^2)
mu.angle
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
                           mu.angle, kappa.angle){
  
  n <- dim(x)[1]
  allprobs <- matrix(1,n,2)
  ind <- which(!is.na(x$angle))[-1] # change: we omit first step 
  # in order to always have the step in t-1
  
  allprobs[ind,] <- cbind(
    dvm(x$angle[ind],mu.angle[1],kappa.angle[1]),
    dvm(x$angle[ind],mu.angle[2],kappa.angle[2])
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
                                 Gamma, delta, autocor, mu.angle, kappa.angle)
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
plot(schwalbe_77$angle, type='l')
points(schwalbe_77$angle, pch=19, col=states_global+1) # vernünftige viz fehlt

hist(schwalbe_77$angle, prob=T, breaks=250, xlab="Step size",
     ylim=c(0,17),xlim=c(-0.5,0.5))

curve(delta[1]*dvm(x,mu.angle[1],kappa.angle[1]), -pi,pi, add=T,n=10000,
      col=4,lwd=2)
curve(delta[2]*dvm(x,mu.angle[2],kappa.angle[2]), -pi,pi, add=T,n=10000,
      col=7,lwd=2)
curve(delta[1]*dvm(x,mu.angle[1],kappa.angle[1])+
        delta[2]*dvm(x,mu.angle[2],kappa.angle[2]), -pi,pi, add=T,n=10000,
      col=2,lwd=2)

# Interessant...





###
### 3-state HMM
###
# starting values
theta = c(
  rep(-2,6), # values to construct TPM
  0.1,0.2,0.2, # autocorrelation [0,1]
  0.1,0.1,0.1, # means of angle for each state [-pi,pi]
  3, 10,15 # kappa of step for each state [0,Inf)
) 
# transformation for unconstrained optimization
theta.star = c(
  theta[1:6], # values to construct TPM
  qlogis(theta[7:9]), # autocorrelation
  theta[10:12] * cos(theta[13:15]), # angle mean
  theta[10:12] * sin(theta[13:15]) # angle kappa
)

mllk_ar1_stepsize(theta.star=theta.star, x=schwalbe_77, 3) # works

# minimize -logL
mod <- nlm(mllk_ar1_stepsize,theta.star,x=schwalbe_77[2:dim(schwalbe_77)[1],],N=3,print.level=2,
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

# angle
mu.angle <- Arg(mod$estimate[(N-1)*N+(N+1):(2*N)]+1i*mod$estimate[(N-1)*N+2*N+1:N])
kappa.angle <- sqrt(mod$estimate[(N-1)*N+(N+1):(2*N)]^2+mod$estimate[(N-1)*N+2*N+1:N]^2)
mu.angle
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
                           mu.angle, kappa.angle){
  
  n <- dim(x)[1]
  allprobs <- matrix(1,n,3)
  ind <- which(!is.na(x$angle))[-c(1:1)] # change: we omit first step 
  # in order to always have the step in t-1
  
  allprobs[ind,] <- cbind(
    dvm(x$angle[ind],mu.angle[1],kappa.angle[1]),
    dvm(x$angle[ind],mu.angle[2],kappa.angle[2]),
    dvm(x$angle[ind],mu.angle[3],kappa.angle[3])
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
                                 Gamma, delta, autocor, mu.angle, kappa.angle)
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
plot(schwalbe_77$angle, type='l')
points(schwalbe_77$angle, pch=19, col=states_global+1) # vernünftige viz fehlt

hist(schwalbe_77$angle, prob=T, breaks=80, xlab="Step size",
     ylim = c(0,12),xlim=c(-pi,pi))

x <- seq(-pi,pi,by=0.0005)
curve(delta[1]*dvm(x,mu.angle[1],kappa.angle[1]), -pi,pi, add=T,
      col=4,lwd=2)
curve(delta[2]*dvm(x,mu.angle[2],kappa.angle[2]), -pi,pi, add=T,
      col=6,lwd=2)
curve(delta[3]*dvm(x,mu.angle[3],kappa.angle[3]), -pi,pi, add=T,
      col=7,lwd=2)
curve(delta[1]*dvm(x,mu.angle[1],kappa.angle[1])+
        delta[2]*dvm(x,mu.angle[2],kappa.angle[2])+
        delta[3]*dvm(x,mu.angle[3],kappa.angle[3]), -pi,pi, add=T,
      col=2,lwd=2)


# ach kacke




###
### next steps: gucken, warum, die einzelnen States so komisch sind...
###





### check for more starting values, takes some time 
# (may need to restart due to singularity issues for some starting values)
llks <- rep(NA,10)
mods <- list()

N=3

for (iteration in 10:10){
  # starting values
  theta = c(
    rep(-2,N*(N-1)), # values to construct TPM
    runif(N,0,1), # autocorrelation [0,1]
    runif(N,-2,2), # means of angle for each state [-pi,pi]
    runif(N,5,1000) # kappa of angle for each state [0,Inf)
  ) 
  # transformation for unconstrained optimization
  theta.star = c(
    theta[1:(N*(N-1))], # values to construct TPM
    qlogis(theta[N*(N-1)+1:N]), # autocorrelation
    # for parameter of von Mises distribution: same transformation as in 
    # function n2w of moveHMM package
    theta[N*(N-1)+N+1:N] * cos(theta[N*(N-1)+2*N+1:N]), # angle mean
    theta[N*(N-1)+N+1:N] * sin(theta[N*(N-1)+2*N+1:N]) # angle kappa
  )
  
  # minimize -logL
  mod <- nlm(mllk_ar1_turningangle,theta.star,x=schwalbe_77[2:4480,],N=2,print.level=0,
             iterlim = 1000)
  mods <- append(mods, list(mod$estimate))
  llks[iteration] <- -mod$minimum
}

llks
mods

sum(round(llks,3)==round(max(llks),3))/length(llks) # proportion of optimal runs
max(llks)

mod$estimate <- mods[[4]]



