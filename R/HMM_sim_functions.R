# 2022-05-30
# Functions for simulating and evaluating HMMs with or without autocorrelation
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # RStudio way


#' Simulate data from a Gamma HMM
#' 
#' @param n_samples Number of samples to generate.
#' @param delta Initial distribution of the Markov chain.
#' @param Gamma Transition probability matrix of the Markov chain.
#' @param N Number of states.
#' @param mu Parameter vector for mu of the gamma distribution.
#' @param sigma Parameter vector for sigma of the gamma distribution.
#' 
#' @return Tupel of states and data of the HMM.
#' 
#' @export
#' @rdname sample_hmm_normal
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


#' Simulate data from a Gamma HMM with AR(1) structure
#' 
#' @param n_samples Number of samples to generate.
#' @param delta Initial distribution of the Markov chain.
#' @param Gamma Transition probability matrix of the Markov chain.
#' @param N Number of states.
#' @param mu Parameter vector for mu of the gamma distribution.
#' @param sigma Parameter vector for sigma of the gamma distribution.
#' @param autocor Parameter vector for the autocorrelation coefficient.
#' 
#' @return Tupel of states and data of the HMM.
#' 
#' @export
#' @rdname sample_hmm_ar1
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


#' Simulate data from a Gamma HMM with AR(p) structure
#' 
#' @param n_samples Number of samples to generate.
#' @param delta Initial distribution of the Markov chain.
#' @param Gamma Transition probability matrix of the Markov chain.
#' @param N Number of states.
#' @param mu Parameter vector for mu of the gamma distribution.
#' @param sigma Parameter vector for sigma of the gamma distribution.
#' @param autocor Parameter vector for the autocorrelation coefficient. 
#'             Has to match p, in the order \eqn{\phi_{t-p},\dots,\phi_{t-1}}
#'             where \eqn{\phi} is the vector of autocorrelation coefficients
#'             for one specific time lag (one value for each state).
#' @param p Degree of autocorrelation.
#' 
#' @return Tupel of states and data of the HMM.
#' 
#' @export
#' @rdname sample_hmm_arp
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


#' Calculate negative Log-Likelihood of gamma HMM
#' 
#' Calculate the negative Log-Likelihood of a \emph{normal} gamma HMM 
#' (without autocorrelation). The forward method using standardized forward 
#' variables is used here.
#' 
#' @param theta.star Unconstrained parameter vector (has to be provided in suitable form).
#' @param x Data vector for which the negative Log-Likelihood should be computed.
#' @param N Number of states.
#' 
#' @return Minus Log-Likelihood
#' 
#' @export
#' @rdname mllk_hmm
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

#' Calculate negative Log-Likelihood of gamma HMM
#' 
#' Calculate the negative Log-Likelihood of an \emph{AR(p)}-gamma HMM 
#' The forward method using standardized forward variables is used here.
#' 
#' @param theta.star Unconstrained parameter vector (has to be provided in suitable form).
#' @param x Data vector for which the negative Log-Likelihood should be computed.
#' @param N Number of states.
#' @param p Degree of autocorrelation.
#' 
#' @return Minus Log-Likelihood
#' 
#' @export
#' @rdname mllk_arp
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


#' Fit an AR(p)-gamma HMM to data
#' 
#' Fit an AR(p)-gamma HMM to data using a specified function to compute the negative 
#' Log-Likelihood. This function gets minimized using the function \code{optim}.
#' It returns the estimated parameters of the fitted model.
#' 
#' @param mllk Negative log-Likelihood function that should be minimized.
#' @param data Data that should be fitted using the HMM.
#' @param theta.star Unconstrained parameter vector (has to be provided in suitable form).
#' @param N Number of states.
#' @param p Degree of autocorrelation, 0 is possible for no autocorrelation.
#' 
#' @return List, containing Gamma, delta, (autocorrelation, depending on degree),
#'         mu, sigma.
#' 
#' @export
#' @rdname fit_arp_model
fit_arp_model <- function(mllk, data, theta.star, N, p){
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

