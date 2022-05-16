##
##
## we use the data from Lawler et al. to test and validate the model 
## formulations
##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set wd


## for the data and their preprocessing we use the markmodmover package 
## from Lawler
library(markmodmover)

sealdata<- data4M(greyseal)
sealdata

sealdata<- interpolate(sealdata,Time.Step = 1.5)
sealdata

seal4M2<- fit(sealdata)
seal4M2


data <- sealdata$data$Movement.Data
head(data)

# data modification to work with our functions
data$Deflection.Angle <- data$Deflection.Angle * pi/180 # degree -> radians
names(data) <- c('angle', 'step', 'group')


#### Model for step size and turning angle with AR(1)-dependence for step size 
###
### 1st aim: Clone model from Lawler et. al (2019)
###       (with exception of von Mises distribution instead of Wrapped Cauchy)
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
  for (t in 2:dim(x)[1]){
    foo <- phi%*%Gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}


# starting values
theta = c(
  -2,-2, # values to construct TPM
  0.2,0.5, # autocorrelation [0,1]
  0.5,1.3, # means of step for each state [0,Inf)
  1, 2, # sd of step for each state [0,Inf)
  0, 0, # means of angle for each state [-pi,pi]
  2, 4 # kappa of angle for each state [0,Inf)
) 
# transformation for unconstrained optimization
theta.star = c(
  theta[1:2], # values to construct TPM
  qlogis(theta[3:4]), # autocorrelation
  log(theta[5:6]), # step mean
  log(theta[7:8]), # step sd
  # for parameter of von Mises distribution: same transformation as in 
  # function n2w of moveHMM package
  theta[9:10] * cos(theta[11:12]), # angle mean
  theta[9:10] * sin(theta[11:12]) # angle kappa
)

mllk_ar1(theta.star=theta.star, x=data, 2) # works

# minimize -logL
mod <- nlm(mllk_ar1,theta.star,x=data[2:dim(data)[1],],N=2,print.level=2,
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
      dvm(x$angle[ind],mu.angle[1],kappa.angle[1]), 
    dgamma(x$step[ind],
           shape=mu.step[2]^2/sigma.step[2]^2,
           scale=sigma.step[2]^2/mu.step[2]) *
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

states_global <- autocor.viterbi(data, 
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
plot(data$step, type='l')
points(data$step, pch=19, col=states_global+1) # vernünftige viz fehlt

hist(data$step, prob=T, breaks=40, xlab="Step size",ylim=c(0,1))

curve(delta[1]*dgamma(x, shape=mu.step[1]^2/sigma.step[1]^2,
                      scale=sigma.step[1]^2/mu.step[1]), 0,45, add=T,n=10000,
      col=5,lwd=2)
curve(delta[2]*dgamma(x, shape=mu.step[2]^2/sigma.step[2]^2,
                      scale=sigma.step[2]^2/mu.step[2]), 0,45,add=T,n=10000,
      col=6,lwd=2)
curve(delta[1]*dgamma(x, shape=mu.step[1]^2/sigma.step[1]^2,
                      scale=sigma.step[1]^2/mu.step[1])+
        delta[2]*dgamma(x, shape=mu.step[2]^2/sigma.step[2]^2,
                      scale=sigma.step[2]^2/mu.step[2]), 0,45,add=T,n=10000,
      col=2,lwd=2)
legend('topright',c('state 1','state 2'),col = c(5,6),pch=15,bty='n')


# angle
plot(data$angle, type='l')
points(data$angle, pch=19, col=states_global+1) # vernünftige viz fehlt

hist(data$angle, prob=T, breaks=70, xlab="angle",xlim=c(-pi,pi))

curve(delta[1]*dvm(x, mu.angle[1],kappa.angle[1]), -pi,pi,add=T,n=10000,
      col=5,lwd=2)
curve(delta[2]*dvm(x ,mu.angle[2],kappa.angle[2]), -pi,pi,add=T,n=10000,
      col=6,lwd=2)
curve(delta[1]*dvm(x, mu.angle[1],kappa.angle[1])+
        delta[2]*dvm(x ,mu.angle[2],kappa.angle[2]), -pi,pi,add=T,n=10000,
      col=2,lwd=2)
legend('topright',c('state 1','state 2'),col = c(5,6),pch=15,bty='n')






### check for more starting values, takes some time 
# (may need to restart due to singularity issues for some starting values)
llks <- rep(NA,100)
mods <- list()

N=2

for (iteration in 1:100){
  # starting values
  theta = c(
    rep(-2,N), # values to construct TPM
    runif(N,0,1), # autocorrelation [0,1]
    runif(N,0.3,1.8), # means of step for each state [0,Inf)
    runif(N,0.1,10), # sd of step for each state [0,Inf)
    runif(N,-2,2), # means of angle for each state [-pi,pi]
    runif(N,5,100) # kappa of angle for each state [0,Inf)
  ) 
  # transformation for unconstrained optimization
  theta.star = c(
    theta[1:2], # values to construct TPM
    qlogis(theta[3:4]), # autocorrelation
    log(theta[5:6]), # step mean
    log(theta[7:8]), # step sd
    # for parameter of von Mises distribution: same transformation as in 
    # function n2w of moveHMM package
    theta[9:10] * cos(theta[11:12]), # angle mean
    theta[9:10] * sin(theta[11:12]) # angle kappa
  )
  
  # minimize -logL
  mod <- nlm(mllk_ar1,theta.star,x=data[2:dim(data)[1],],N=2,print.level=0,
             iterlim = 1000)
  append(mods, mod)
  llks[iteration] <- -mod$minimum
}

llks

sum(round(llks,3)==round(max(llks),3))/length(llks) # proportion of optimal runs
max(llks)


#####
##### N=3
#####
#####

# starting values
N=3
theta = c(
  rep(-2,N*(N-1)), # values to construct TPM
  runif(N,0,1), # autocorrelation [0,1]
  runif(N,0.2,2.5), # means of step for each state [0,Inf)
  runif(N,0.1,10), # sd of step for each state [0,Inf)
  runif(N,-2,2), # means of angle for each state [-pi,pi]
  runif(N,5,100) # kappa of angle for each state [0,Inf)
) 
# transformation for unconstrained optimization
theta.star = c(
  theta[1:(N*(N-1))], # values to construct TPM
  qlogis(theta[N*(N-1)+1:N]), # autocorrelation
  log(theta[N*(N-1)+N+1:N]), # step mean
  log(theta[N*(N-1)+2*N+1:N]), # step sd
  # for parameter of von Mises distribution: same transformation as in 
  # function n2w of moveHMM package
  theta[N*(N-1)+3*N+1:N] * cos(theta[N*(N-1)+4*N+1:N]), # angle mean
  theta[N*(N-1)+3*N+1:N] * sin(theta[N*(N-1)+4*N+1:N]) # angle kappa
)

# minimize -logL
mod <- nlm(mllk_ar1,theta.star,x=data[2:dim(data)[1],],N=3,print.level=2,
           iterlim = 1000)



## re-transformation to natural parameters
## TPM
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
### state decoding
###
###
### global decoding, using viterbi from lecture
###

# We need TPM Gamma
autocor.viterbi <-function(x, Gamma, delta, autocor, 
                           mu.step, sigma.step, mu.angle, kappa.angle){
  
  n <- dim(x)[1]
  allprobs <- matrix(1,n,3)
  ind <- which(!is.na(x$step))[-1] # change: we omit first step 
  # in order to always have the step in t-1
  
  allprobs[ind,] <- cbind(
    dgamma(x$step[ind],
           shape=mu.step[1]^2/sigma.step[1]^2,
           scale=sigma.step[1]^2/mu.step[1]) *
      dvm(x$angle[ind],mu.angle[1],kappa.angle[1]), 
    dgamma(x$step[ind],
           shape=mu.step[2]^2/sigma.step[2]^2,
           scale=sigma.step[2]^2/mu.step[2]) *
      dvm(x$angle[ind],mu.angle[2],kappa.angle[2]),
    dgamma(x$step[ind],
           shape=mu.step[3]^2/sigma.step[3]^2,
           scale=sigma.step[3]^2/mu.step[3]) *
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

states_global <- autocor.viterbi(data, 
                                 Gamma, delta, autocor, mu.step, sigma.step,
                                 mu.angle, kappa.angle)
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
plot(data$step, type='l')
points(data$step, pch=19, col=states_global+1) # vernünftige viz fehlt

hist(data$step, prob=T, breaks=40, xlab="Step size",ylim=c(0,1))

curve(delta[1]*dgamma(x, shape=mu.step[1]^2/sigma.step[1]^2,
                      scale=sigma.step[1]^2/mu.step[1]), 0,45, add=T,n=10000,
      col=5,lwd=2)
curve(delta[2]*dgamma(x, shape=mu.step[2]^2/sigma.step[2]^2,
                      scale=sigma.step[2]^2/mu.step[2]), 0,45,add=T,n=10000,
      col=6,lwd=2)
curve(delta[3]*dgamma(x, shape=mu.step[3]^2/sigma.step[3]^2,
                      scale=sigma.step[3]^2/mu.step[3]), 0,45,add=T,n=10000,
      col=7,lwd=2)
curve(delta[1]*dgamma(x, shape=mu.step[1]^2/sigma.step[1]^2,
                      scale=sigma.step[1]^2/mu.step[1])+
        delta[2]*dgamma(x, shape=mu.step[2]^2/sigma.step[2]^2,
                        scale=sigma.step[2]^2/mu.step[2])+
        delta[3]*dgamma(x, shape=mu.step[3]^2/sigma.step[3]^2,
                        scale=sigma.step[3]^2/mu.step[3]), 0,45,add=T,n=10000,
      col=2,lwd=2)
legend('topright',c('state 1','state 2', 'state 3'),col = c(5,6,7),pch=15,bty='n')


# angle
plot(data$angle, type='l')
points(data$angle, pch=19, col=states_global+1) # vernünftige viz fehlt

hist(data$angle, prob=T, breaks=70, xlab="angle",xlim=c(-pi,pi))

curve(delta[1]*dvm(x, mu.angle[1],kappa.angle[1]), -pi,pi,add=T,n=10000,
      col=5,lwd=2)
curve(delta[2]*dvm(x ,mu.angle[2],kappa.angle[2]), -pi,pi,add=T,n=10000,
      col=6,lwd=2)
curve(delta[3]*dvm(x ,mu.angle[3],kappa.angle[3]), -pi,pi,add=T,n=10000,
      col=7,lwd=2)
curve(delta[1]*dvm(x, mu.angle[1],kappa.angle[1])+
        delta[2]*dvm(x ,mu.angle[2],kappa.angle[2])+
        delta[3]*dvm(x ,mu.angle[3],kappa.angle[3]), -pi,pi,add=T,n=10000,
      col=2,lwd=2)
legend('topright',c('state 1','state 2','state 3'),col = c(5,6,7),pch=15,bty='n')







