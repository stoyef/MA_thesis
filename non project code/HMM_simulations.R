## 2022-05-23
## HMM-Simulations


#### Goal: Simulate from different kinds of HMMs and compare the performance
####       by refitting 
### - Normal HMM
### - HMM with AR(1)-process in state dependent distribution
### - HMM with higher order AR-processes in state dependent distribution

### We need to specify: Kind of distribution for states, number of states,
###                     Initial distribution, TPM

##
## (Initial) assumption: univariate Gamma distribution, 2 states
##


################################################################################
##################### Normal HMM ###############################################
################################################################################

sample_hmm_normal <- function(n_samples, delta, Gamma, N, mu, sigma){
  states <- rep(NA,n_samples)
  data <- rep(NA,n_samples)
  
  states[1] <- sample(1:N, 1, prob = delta)
  data[1] <- rgamma(1,shape=mu[states[1]]^2/sigma[states[1]]^2,
                    scale=sigma[states[1]]^2/mu[states[1]])
  
  for (t in 2:n_samples){
    states[t] <- sample(1:N, 1, prob=Gamma[states[t-1],])
    data[t] <- rgamma(1,shape=mu[states[t]]^2/sigma[states[t]]^2,
                      scale=sigma[states[t]]^2/mu[states[t]])
  }
  
  ret <- list(states, data)
  names(ret) <- c('states','data')
  return(ret)
}

delta_normal <- c(0.5,0.5)
Gamma_normal <- matrix(c(0.8,0.2,0.2,0.8),nrow=2)
mu_normal <- c(10,20)
sigma_normal <- c(2,5)
sim_normal <- sample_hmm_normal(1000, delta_normal, Gamma_normal, 2, 
                                mu_normal, sigma_normal)
sim_normal$data


################################################################################
####################### HMM with AR(1) #########################################
################################################################################

# Attention: We can only compute the HMM given the first data point (due to AR(1))

sample_hmm_ar1 <- function(n_samples, delta, Gamma, N, mu, sigma, autocor){
  states <- rep(NA,n_samples)
  data <- rep(NA,n_samples)
  
  # given, no autocorrelation -> just like normal HMM
  states[1] <- sample(1:N, 1, prob = delta)
  data[1] <- rgamma(1,shape=mu[states[1]]^2/sigma[states[1]]^2,
                    scale=sigma[states[1]]^2/mu[states[1]])
  
  for (t in 2:n_samples){
    states[t] <- sample(1:N, 1, prob=Gamma[states[t-1],])
    # compute mu with AR(1)
    mu_ar <- autocor[states[t]]*data[t-1] + (1-autocor[states[t]])*mu[states[t]]
    data[t] <- rgamma(1,shape=mu_ar^2/sigma[states[t]]^2,
                  scale=sigma[states[t]]^2/mu_ar)
  }
  
  ret <- list(states, data)
  names(ret) <- c('states','data')
  return(ret)
}

delta_ar1 <- c(0.5,0.5)
Gamma_ar1 <- matrix(c(0.8,0.2,0.2,0.8),nrow=2)
mu_ar1 <- c(10,20)
sigma_ar1 <- c(2,5)
autocor_ar1 <- c(0.2,0.4)
sim_ar1 <- sample_hmm_ar1(1000, delta_ar1, Gamma_ar1, 2, 
                                mu_ar1, sigma_ar1, autocor_ar1)
sim_ar1$data


################################################################################
####################### HMM with AR(p) #########################################
################################################################################

# Attention: We can only compute the HMM given the first p data points (due to AR(p))
# autocor has to be supplied as matrix of the autocorrelation coefficients, with
# states in rows and time lags in columns!

sample_hmm_arp <- function(n_samples, delta, Gamma, N, mu, sigma, autocor, p){
  states <- rep(NA,n_samples)
  data <- rep(NA,n_samples)
  
  # first p data points given, no autocorrelation -> just like normal HMM
  states[1] <- sample(1:N, 1, prob = delta)
  data[1] <- rgamma(1,shape=mu[states[1]]^2/sigma[states[1]]^2,
                    scale=sigma[states[1]]^2/mu[states[1]])
  if (p>1){
    for (t in 2:p){
    states[t] <- sample(1:N, 1, prob = Gamma[states[t-1],])
    data[t] <- rgamma(1,shape=mu[states[t]]^2/sigma[states[t]]^2,
                      scale=sigma[states[t]]^2/mu[states[t]])
    }
  }
  
  for (t in (p+1):n_samples){
    states[t] <- sample(1:N, 1, prob=Gamma[states[t-1],])
    # compute mu with AR(p)
    ##
    ###
    ## hier: muss summe der autocor 1 ergeben oder alles in der Summe 1????
    ## -> wir machen erstmal alles in der Summe 1 aber nochmal überlegen!
    mu_ar <- sum(autocor[states[t],]*data[(t-p):(t-1)]) + 
      (1-sum(autocor[states[t],]))*mu[states[t]]
    data[t] <- rgamma(1,shape=mu_ar^2/sigma[states[t]]^2,
                      scale=sigma[states[t]]^2/mu_ar)
  }
  
  ret <- list(states, data)
  names(ret) <- c('states','data')
  return(ret)
}

delta_arp <- c(0.5,0.5)
Gamma_arp <- matrix(c(0.8,0.2,0.2,0.8),nrow=2)
mu_arp <- c(10,20)
sigma_arp <- c(2,5)
autocor_arp <- matrix(c(0.1,0.2,0.3,0.4),nrow=2) # first entry for t-2, second entry for t-1
p <- 2
sim_arp <- sample_hmm_arp(1000, delta_arp, Gamma_arp, 2, 
                          mu_arp, sigma_arp, autocor_arp, p)
sim_arp$data


################################################################################
################################################################################
################################################################################
# seems to work
# now we apply each model to each simulated data frame


# minus log-likelihood normal HMM
mllk_hmm <-function(theta.star,x,N){
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
                                  shape=mu[j]^2/sigma[j]^2,
                                  scale=sigma[j]^2/mu[j]) 
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


# minus log-likelihood AR(p) HMM
mllk_arp <-function(theta.star,x,N,p){
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  
  autocor <- plogis(theta.star[(N-1)*N+1:(p*N)])
  autocor <- matrix(autocor, ncol=p, byrow=TRUE) # matrix for easier handling later on
  mu <- exp(theta.star[(N-1)*N+(p*N+1):(p*N+N)])
  sigma <- exp(theta.star[(N-1)*N+p*N+N+1:N])
  
  allprobs <- matrix(1,length(x),N)
  ind <- which(!is.na(x))[-c(1:p)] # change: we omit first p steps 
  # in order to always have the step in t-p
  
  autocor_ind <- matrix(NA,nrow=length(ind),ncol=p) # matrix for indices of autocor data
  for (i in 1:p){
    autocor_ind[,i] <- ind-p+i-1
  }
  autocor_ind <- apply(autocor_ind, 2, function(a)x[a]) # substitute indices with values
  
  for (j in 1:N){
    # here comes the autocorrelation!
    mu_auto <- c(rep(NA,p), # AR(p)
                 ((1-sum(autocor[j,]))*mu[j] + 
                    as.vector(autocor_ind%*%autocor[j,]))) # matmul of values with autocor coefficient
    allprobs[ind,j] <- dgamma(x[ind],
                                  shape=mu_auto[ind]^2/sigma[j]^2,
                                  scale=sigma[j]^2/mu_auto[ind]) 
    # here we have to choose mu_auto[ind], because
    # we have an individual mu for each data point
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

# Again, the autocorrelation parameters are supplied as a matrix. To make handling
# easier, we supply a vector in the input which then gets converted into a 
# matrix. Attention: The matrix is filled row-wise in reverse order of the time lags
# Meaning: First p values in vector become the first row, a[1,1] corresponds to 
# t-p of state 1, a[1,2] corresponds to t-p+1 of state 1 and so on
# -> Matrix has dimensions Nxp

# For all data sets we choose as many starting parameters the same as possible

theta_normal <- c(
  -2,-2, # TPM
  15,30, # mu
  5,10 # sigma
)

theta_normal_star <- c(
  theta_normal[1:2],
  log(theta_normal[3:6])
)

theta_ar1 <- c(
  -2,-2, # TPM
  0.2,0.4, # autocorrelation
  15,30, # mu
  5,10 # sigma
)

theta_ar1_star <- c(
  theta_ar1[1:2],
  qlogis(theta_ar1[3:4]),
  log(theta_ar1[5:8])
)

theta_arp <- c(
  -2,-2, # TPM
  0.1,0.1,0.3,0.1, # autocorrelation, in this case data was generated from AR(2)
  15,30, # mu
  5,10 # sigma
)

theta_arp_star <- c(
  theta_arp[1:2],
  qlogis(theta_arp[3:6]),
  log(theta_arp[7:10])
)

###
### Data from normal HMM
###

# Minimize -logL

####### Normal HMM #######
mod_hmm <- nlm(mllk_hmm,theta_normal_star,x=sim_normal$data,N=2,print.level=2,
           iterlim=1000)
mod_hmm$estimate
## re-transformation to natural parameters
## TPM
N=2
Gamma <- diag(N)
Gamma[!Gamma] <- exp(mod_hmm$estimate[1:((N-1)*N)])
Gamma <- Gamma/rowSums(Gamma)
Gamma
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
delta
# mu and sigma 
mu <- exp(mod_hmm$estimate[(N-1)*N+1:N])
mu
sigma <- exp(mod_hmm$estimate[(N-1)*N+(N+1):(2*N)])
sigma


####### AR(1)-HMM ########
mllk_arp(theta_ar1_star,sim_normal$data,N=2,p=1)
mod_ar1 <- nlm(mllk_arp,theta_ar1_star,x=sim_normal$data,N=2,p=1,print.level=2,
           iterlim=1000)
## nlm doesn't work here for some reason, let's try optim with L-BFGS instead
mod_ar1 <- optim(par=theta_ar1_star, fn=mllk_arp, method='L-BFGS-B',
                 N=2,p=1,x=sim_normal$data)
mod_ar1$par
# yay, seems to work

## re-transformation to natural parameters
## TPM
N=2
Gamma <- diag(N)
Gamma[!Gamma] <- exp(mod_ar1$par[1:((N-1)*N)])
Gamma <- Gamma/rowSums(Gamma)
Gamma
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
delta
# autocorrelation
plogis(mod_ar1$par[(N-1)*N+1:N])
# mu and sigma 
mu <- exp(mod_ar1$par[(N-1)*N+(N+1):(2*N)])
mu
sigma <- exp(mod_ar1$par[(N-1)*N+2*N+1:N])
sigma
# looking good :)


####### AR(2)-HMM ########
mllk_arp(theta_arp_star,sim_normal$data,N=2,p=2)
mod_arp <- optim(par=theta_arp_star, fn=mllk_arp, method='L-BFGS-B',
                 N=2,p=2,x=sim_normal$data)
mod_arp$par
## re-transformation to natural parameters
## TPM
N=2
p=2
Gamma <- diag(N)
Gamma[!Gamma] <- exp(mod_arp$par[1:((N-1)*N)])
Gamma <- Gamma/rowSums(Gamma)
Gamma
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
delta
# autocorrelation
plogis(mod_arp$par[(N-1)*N+1:(p*N)])
# mu and sigma 
mu <- exp(mod_arp$par[(N-1)*N+(p*N+1):(p*N+N)])
mu
sigma <- exp(mod_arp$par[(N-1)*N+p*N+N+1:N])
sigma




###
### Data from AR(1) HMM
###

# Minimize -logL

####### Normal HMM #######
mod_hmm <- nlm(mllk_hmm,theta_normal_star,x=sim_ar1$data,N=2,print.level=2,
               iterlim=1000)
mod_hmm$estimate
## re-transformation to natural parameters
## TPM
N=2
Gamma <- diag(N)
Gamma[!Gamma] <- exp(mod_hmm$estimate[1:((N-1)*N)])
Gamma <- Gamma/rowSums(Gamma)
Gamma
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
delta
# mu and sigma 
mu <- exp(mod_hmm$estimate[(N-1)*N+1:N])
mu
sigma <- exp(mod_hmm$estimate[(N-1)*N+(N+1):(2*N)])
sigma


####### AR(1)-HMM ########
mllk_arp(theta_ar1_star,sim_ar1$data,N=2,p=1)
mod_ar1 <- optim(par=theta_ar1_star, fn=mllk_arp, method='L-BFGS-B',
                 N=2,p=1,x=sim_ar1$data)
mod_ar1$par
## re-transformation to natural parameters
## TPM
N=2
Gamma <- diag(N)
Gamma[!Gamma] <- exp(mod_ar1$par[1:((N-1)*N)])
Gamma <- Gamma/rowSums(Gamma)
Gamma
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
delta
# autocorrelation
plogis(mod_ar1$par[(N-1)*N+1:N])
# mu and sigma 
mu <- exp(mod_ar1$par[(N-1)*N+(N+1):(2*N)])
mu
sigma <- exp(mod_ar1$par[(N-1)*N+2*N+1:N])
sigma


####### AR(2)-HMM ########
mllk_arp(theta_arp_star,sim_arp$data,N=2,p=2)
mod_arp <- optim(par=theta_arp_star, fn=mllk_arp, method='L-BFGS-B',
                 N=2,p=2,x=sim_ar1$data)
mod_arp$par
## re-transformation to natural parameters
## TPM
N=2
p=2
Gamma <- diag(N)
Gamma[!Gamma] <- exp(mod_arp$par[1:((N-1)*N)])
Gamma <- Gamma/rowSums(Gamma)
Gamma
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
delta
# autocorrelation
plogis(mod_arp$par[(N-1)*N+1:(p*N)])
# mu and sigma 
mu <- exp(mod_arp$par[(N-1)*N+(p*N+1):(p*N+N)])
mu
sigma <- exp(mod_arp$par[(N-1)*N+p*N+N+1:N])
sigma




###
### Data from AR(2) HMM
###

# Minimize -logL

####### Normal HMM #######
mod_hmm <- nlm(mllk_hmm,theta_normal_star,x=sim_arp$data,N=2,print.level=2,
               iterlim=1000)
mod_hmm$estimate
## re-transformation to natural parameters
## TPM
N=2
Gamma <- diag(N)
Gamma[!Gamma] <- exp(mod_hmm$estimate[1:((N-1)*N)])
Gamma <- Gamma/rowSums(Gamma)
Gamma
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
delta
# mu and sigma 
mu <- exp(mod_hmm$estimate[(N-1)*N+1:N])
mu
sigma <- exp(mod_hmm$estimate[(N-1)*N+(N+1):(2*N)])
sigma


####### AR(1)-HMM ########
mllk_arp(theta_ar1_star,sim_ar1$data,N=2,p=1)
mod_ar1 <- optim(par=theta_ar1_star, fn=mllk_arp, method='L-BFGS-B',
                 N=2,p=1,x=sim_arp$data)
mod_ar1$par
## re-transformation to natural parameters
## TPM
N=2
Gamma <- diag(N)
Gamma[!Gamma] <- exp(mod_ar1$par[1:((N-1)*N)])
Gamma <- Gamma/rowSums(Gamma)
Gamma
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
delta
# autocorrelation
plogis(mod_ar1$par[(N-1)*N+1:N])
# mu and sigma 
mu <- exp(mod_ar1$par[(N-1)*N+(N+1):(2*N)])
mu
sigma <- exp(mod_ar1$par[(N-1)*N+2*N+1:N])
sigma


####### AR(2)-HMM ########
mllk_arp(theta_arp_star,sim_arp$data,N=2,p=2)
mod_arp <- optim(par=theta_arp_star, fn=mllk_arp, method='L-BFGS-B',
                 N=2,p=2,x=sim_arp$data)
mod_arp$par
## re-transformation to natural parameters
## TPM
N=2
p=2
Gamma <- diag(N)
Gamma[!Gamma] <- exp(mod_arp$par[1:((N-1)*N)])
Gamma <- Gamma/rowSums(Gamma)
Gamma
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
delta
# autocorrelation
plogis(mod_arp$par[(N-1)*N+1:(p*N)])
# mu and sigma 
mu <- exp(mod_arp$par[(N-1)*N+(p*N+1):(p*N+N)])
mu
sigma <- exp(mod_arp$par[(N-1)*N+p*N+N+1:N])
sigma


## quite good results, next step: Quantitative evaluation and global decoding


## for easier use: function to fit model

fit.model <- function(mllk, data, theta.star, N, p){
  if (p>0){
  mod <- optim(par=theta.star, fn=mllk, method='L-BFGS-B',
               N=N,p=p,x=data)
  } else{
  mod <- optim(par=theta.star, fn=mllk, method='L-BFGS-B',
                 N=N,x=data)
  }
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(mod$par[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  if (p>0){
  # autocorrelation
  autocor <- plogis(mod$par[(N-1)*N+1:(p*N)])
  # mu and sigma 
  mu <- exp(mod$par[(N-1)*N+(p*N+1):(p*N+N)])
  sigma <- exp(mod$par[(N-1)*N+p*N+N+1:N])
  ret <- list(Gamma, delta, autocor, mu, sigma)
  names(ret) <- c('Gamma', 'delta', 'autocorrelation', 'mu', 'sigma')
  return(ret)
  } else{
    # mu and sigma 
    mu <- exp(mod$par[(N-1)*N+1:N])
    sigma <- exp(mod$par[(N-1)*N+N+1:N])
    ret <- list(Gamma, delta, mu, sigma)
    names(ret) <- c('Gamma', 'delta', 'mu', 'sigma')
    return(ret)
  }
}



mod_arp <- fit.model(mllk_arp, sim_arp$data, theta_arp_star, N=2,p=2)
mod_arp

mod_ar1 <- fit.model(mllk_arp, sim_arp$data, theta_ar1_star, N=2,p=1)
mod_ar1

mod_hmm <- fit.model(mllk_hmm, sim_arp$data, theta_normal_star, N=2,p=0)
mod_hmm


#######
#######
####### "Model selection" using AIC, BIC
#######
#######

theta_arp_best_arp <- c(
  mod_arp$Gamma[1,2],mod_arp$Gamma[2,1], # TPM
  mod_arp$autocorrelation, # autocorrelation, in this case data was generated from AR(2)
  mod_arp$mu, # mu
  mod_arp$sigma # sigma
)

theta_arp_best_star_arp <- c(
  theta_arp_best[1:2],
  qlogis(theta_arp_best[3:6]),
  log(theta_arp_best[7:10])
)

theta_ar1_best_arp <- c(
  mod_ar1$Gamma[1,2],mod_ar1$Gamma[2,1], # TPM
  mod_ar1$autocorrelation, # autocorrelation, in this case data was generated from AR(2)
  mod_ar1$mu, # mu
  mod_ar1$sigma # sigma
)

theta_ar1_best_star_arp <- c(
  theta_ar1_best[1:2],
  qlogis(theta_ar1_best[3:4]),
  log(theta_ar1_best[5:8])
)

theta_hmm_best_arp <- c(
  mod_hmm$Gamma[1,2],mod_hmm$Gamma[2,1], # TPM
  mod_hmm$mu, # mu
  mod_hmm$sigma # sigma
)

theta_hmm_best_star_arp <- c(
  theta_ar1_best[1:2],
  log(theta_ar1_best[3:6])
)

## AIC, data from AR(2) simulation
# Normal HMM model
2*mllk_hmm(theta_hmm_best_star, sim_arp$data, N=2) + 2*length(theta_hmm_best_star)
# AR(1) model
2*mllk_arp(theta_ar1_best_star, sim_arp$data, N=2,p=1) + 2*length(theta_ar1_best_star)
# AR(2) model
2*mllk_arp(theta_arp_best_star, sim_arp$data, N=2,p=2) + 2*length(theta_arp_best_star)





#######
#######
####### Global decoding
#######
#######


# We decode all fitted models and compare the decoding errors using viterbi


##
##
## Ändert sich bei Viterbi irgendwas???
##
##



viterbi_arp <-function(x, Gamma, delta, autocor, 
                           mu, sigma, N, p){
  n <- length(x)
  allprobs <- matrix(1,n,N)
  ind <- which(!is.na(x))[-c(1:p)] # change: we omit first step 
  # in order to always have the step in t-1
  
  autocor <- matrix(autocor, ncol=p, byrow=TRUE) # aurocorrelation matrix for easier handling later on
  
  autocor_ind <- matrix(NA,nrow=length(ind),ncol=p) # matrix for indices of autocor data
  for (i in 1:p){
    autocor_ind[,i] <- ind-p+i-1
  }
  autocor_ind <- apply(autocor_ind, 2, function(a)x[a]) # substitute indices with values
  
  for (j in 1:N){
  mu_auto <- c(rep(NA,p), # AR(p)
               ((1-sum(autocor[j,]))*mu[j] + 
                  as.vector(autocor_ind%*%autocor[j,]))) # matmul of values with autocor coefficient
  
  allprobs[ind,j] <- dgamma(x[ind],
           shape=mu_auto[ind]^2/sigma[j]^2,
           scale=sigma[j]^2/mu_auto[ind])
  }
  
  xi <- matrix(0,n,N)
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

