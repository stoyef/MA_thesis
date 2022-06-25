## old functions from r-project that are no longer needed


#' Simulate data from a Gamma HMM with AR(p) structure
#' 
#' Simulate data from a Gamma HMM with AR(p) structure in the parameters \eqn{\mu}
#' and \eqn{\sigma}.
#'
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
#'             for one specific time lag (one value for each state),
#'             0 if no autocorrelation.
#' @param p Degree of autocorrelation, 0 = no autocorrelation.
#' 
#' @return Tupel of states and data of the HMM.
#' 
#' @export
#' @rdname sample_gamma_arp
sample_gamma_arp <- function(n_samples, delta, Gamma, N, mu, sigma, autocor, p){
  states <- rep(NA,n_samples)
  data <- rep(NA,n_samples)
  cv <- sigma/mu
  
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
  
  if (p==0){ # no autocorrelation
    for (t in 2:n_samples){
      states[t] <- sample(1:N, 1, prob=Gamma[states[t-1],])
      data[t] <- rgamma(1,shape=mu[states[t]]^2/sigma[states[t]]^2,
                        scale=sigma[states[t]]^2/mu[states[t]])
    }
  } else{
    for (t in (p+1):n_samples){
      states[t] <- sample(1:N, 1, prob=Gamma[states[t-1],])
      # compute mu with AR(p)
      ##
      ###
      ## hier: muss summe der autocor 1 ergeben oder alles in der Summe 1????
      ## -> wir machen erstmal alles in der Summe 1 aber nochmal überlegen!
      mu_ar <- sum(autocor[states[t],]*data[(t-p):(t-1)]) + 
        (1-sum(autocor[states[t],]))*mu[states[t]]
      sigma_ar <- cv[states[t]]*mu_ar # retain ccv
      data[t] <- rgamma(1,shape=mu_ar^2/sigma_ar^2,
                        scale=sigma_ar^2/mu_ar)
    }
  }
  
  ret <- list(states, data)
  names(ret) <- c('states','data')
  return(ret)
}


#' Simulate data from a von Mises HMM with AR(p) structure
#' 
#' Simulate data from a von Mises HMM with AR(p) structure in the parameter \eqn{\mu}
#' (and \eqn{\kappa} in the future?).
#'
#' 
#' @param n_samples Number of samples to generate.
#' @param delta Initial distribution of the Markov chain.
#' @param Gamma Transition probability matrix of the Markov chain.
#' @param N Number of states.
#' @param mu Parameter vector for mu of the gamma distribution.
#' @param kappa Parameter vector for sigma of the gamma distribution.
#' @param autocor Parameter vector for the autocorrelation coefficient. 
#'             Has to match p, in the order \eqn{\phi_{t-p},\dots,\phi_{t-1}}
#'             where \eqn{\phi} is the vector of autocorrelation coefficients
#'             for one specific time lag (one value for each state).
#'             0, if no autocorrelation
#' @param p Degree of autocorrelation, 0=no autocorrelation.
#' 
#' @return Tupel of states and data of the HMM.
#' 
#' @export
#' @rdname sample_vonMises_arp
sample_vonMises_arp <- function(n_samples, delta, Gamma, N, mu, kappa, autocor, p){
  require(CircStats)
  states <- rep(NA,n_samples)
  data <- rep(NA,n_samples)
  
  # first p data points given, no autocorrelation -> just like normal HMM
  # we have to be careful with rvm, it uses 0-2*pi instead of -pi-+pi
  states[1] <- sample(1:N, 1, prob = delta)
  data[1] <- rvm(1,mean=mu[states[1]]+pi,
                 k=kappa[states[1]])-pi
  if (p>1){
    for (t in 2:p){
      states[t] <- sample(1:N, 1, prob = Gamma[states[t-1],])
      data[t] <- rvm(1,mean=mu[states[t]]+pi,
                     k=kappa[states[t]])-pi
    }
  }
  
  if (p==0){ # no autocorrelation
    for (t in 2:n_samples){
      states[t] <- sample(1:N, 1, prob=Gamma[states[t-1],])
      data[t] <- rvm(1,mean=mu[states[t]]+pi,k=kappa[states[t]])-pi
    }
  } else{
    for (t in (p+1):n_samples){
      states[t] <- sample(1:N, 1, prob=Gamma[states[t-1],])
      # compute mu with AR(p)
      ##
      ###
      ## hier: muss summe der autocor 1 ergeben oder alles in der Summe 1????
      ## -> wir machen erstmal alles in der Summe 1 aber nochmal überlegen!
      mu_ar <- sum(autocor[states[t],]*data[(t-p):(t-1)]) + 
        (1-sum(autocor[states[t],]))*mu[states[t]]
      # rvm samples in region 0-2pi, so we substract pi
      data[t] <- rvm(1,mean=mu_ar+pi,
                     k=kappa[states[t]]) - pi
    }
  }
  
  ret <- list(states, data)
  names(ret) <- c('states','data')
  return(ret)
}




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
#' Simulate data from a Gamma HMM with AR(1) structure in the parameter \eqn{\mu}.
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
#' @rdname sample_hmm_ar1_mu
sample_hmm_ar1_mu <- function(n_samples, delta, Gamma, N, mu, sigma, autocor){
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


#' Simulate data from a Gamma HMM with AR(1) structure
#' 
#' Simulate data from a Gamma HMM with AR(1) structure in parameters \eqn{\mu} and \eqn{\sigma}.
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
  cv <- sigma/mu
  
  # given, no autocorrelation -> just like normal HMM
  states[1] <- sample(1:N, 1, prob = delta)
  data[1] <- rgamma(1,shape=mu[states[1]]^2/sigma[states[1]]^2,
                    scale=sigma[states[1]]^2/mu[states[1]])
  
  for (t in 2:n_samples){
    states[t] <- sample(1:N, 1, prob=Gamma[states[t-1],])
    # compute mu with AR(1)
    mu_ar <- autocor[states[t]]*data[t-1] + (1-autocor[states[t]])*mu[states[t]]
    sigma_ar <- cv[states[t]]*mu_ar # retain ccv
    data[t] <- rgamma(1,shape=mu_ar^2/sigma_ar^2,
                      scale=sigma_ar^2/mu_ar)
  }
  
  ret <- list(states, data)
  names(ret) <- c('states','data')
  return(ret)
}


#' Simulate data from a Gamma HMM with AR(p) structure
#' 
#' Simulate data from a Gamma HMM with AR(p) structure in the parameter \eqn{\mu}.

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
#' @rdname sample_hmm_arp_mu
sample_hmm_arp_mu <- function(n_samples, delta, Gamma, N, mu, sigma, autocor, p){
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
#' @param p Degree of autocorrelation, here always = 0.
#' 
#' @return Minus Log-Likelihood
#' 
#' @export
#' @rdname mllk_gamma_hmm
mllk_gamma_hmm <-function(theta.star,x,N,p=0){
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


#' Calculate negative Log-Likelihood of AR(p)-gamma HMM
#' 
#' Calculate the negative Log-Likelihood of an gamma HMM with \emph{AR(p)} only in the
#' parameter \eqn{\mu}.
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
#' @rdname mllk_arp_mu
mllk_arp_mu <-function(theta.star,x,N,p){
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  
  mu <- exp(theta.star[(N-1)*N+1:N])
  sigma <- exp(theta.star[(N-1)*N+N+1:N])
  
  autocor <- plogis(theta.star[(N-1)*N+2*N+1:(p*N)])
  autocor <- matrix(autocor, ncol=p, byrow=TRUE) # matrix for easier handling later on
  
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


#' Calculate negative Log-Likelihood of AR(p)-gamma HMM
#' 
#' Calculate the negative Log-Likelihood of an gamma HMM with \emph{AR(p)} in the
#' parameters \eqn{\mu} and \eqn{\sigma}.
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
#' @rdname mllk_gamma_arp
mllk_gamma_arp <-function(theta.star,x,N,p){
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  
  mu <- exp(theta.star[(N-1)*N+1:N])
  sigma <- exp(theta.star[(N-1)*N+N+1:N])
  # calculate coefficient of variance, should be constant while modeling
  cv <- sigma/mu
  allprobs <- matrix(1,length(x),N)
  
  if (p>0){
    autocor <- plogis(theta.star[(N-1)*N+2*N+1:(p*N)])
    autocor <- matrix(autocor, ncol=p, byrow=TRUE) # matrix for easier handling later on
    
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
      sigma_auto <- cv[j]*mu_auto # calculate sigma using ccv
      allprobs[ind,j] <- dgamma(x[ind],
                                shape=mu_auto[ind]^2/sigma_auto[ind]^2,
                                scale=sigma_auto[ind]^2/mu_auto[ind]) 
      # here we have to choose mu_auto[ind], because
      # we have an individual mu for each data point
    }
  } else{
    ind <- which(!is.na(x))
    
    for (j in 1:N){
      allprobs[ind,j] <- dgamma(x[ind],
                                shape=mu[j]^2/sigma[j]^2,
                                scale=sigma[j]^2/mu[j]) 
    }
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



#' Calculate negative Log-Likelihood of AR(p)-von Mises HMM
#' 
#' Calculate the negative Log-Likelihood of an von Mises HMM with \emph{AR(p)} in the
#' parameter \eqn{\mu} (and \eqn{\kappa} in the future?).
#' The forward method using standardized forward variables is used here.
#' 
#' @param theta.star Unconstrained parameter vector (has to be provided in suitable form:
#' Gamma, mu, kappa, autocorrelation)!
#' @param x Data vector for which the negative Log-Likelihood should be computed.
#' @param N Number of states.
#' @param p Degree of autocorrelation, 0=no autocorrelation.
#' 
#' @return Minus Log-Likelihood
#' 
#' @export
#' @rdname mllk_vonMises_arp
mllk_vonMises_arp <-function(theta.star,x,N,p){
  require(CircStats)
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  
  # same transformation as in function w2n of moveHMM package:
  mu <- Arg(theta.star[(N-1)*N+1:N]+1i*theta.star[(N-1)*N+N+1:N])
  kappa <- sqrt(theta.star[(N-1)*N+1:N]^2+theta.star[(N-1)*N+N+1:N]^2)
  
  allprobs <- matrix(1,length(x),N)
  
  if (p>0){
    autocor <- plogis(theta.star[(N-1)*N+(2*N)+1:(p*N)])
    autocor <- matrix(autocor, ncol=p, byrow=TRUE) # matrix for easier handling later on
    
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
      allprobs[ind,j] <- dvm(x[ind],
                             mu=mu_auto[ind],
                             kappa=kappa[j])
      # here we have to choose mu_auto[ind], because
      # we have an individual mu for each data point
    }
  } else{
    ind <- which(!is.na(x))
    
    for (j in 1:N){
      allprobs[ind,j] <- dvm(x[ind],
                             mu=mu[j],
                             kappa=kappa[j])
    }
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



#' Global decoding for AR(p)-gamma HMM using Viterbi
#' 
#' Global decoding for an AR(p)-gamma HMM with autocorrelation in the parameter \eqn{\mu}
#' using the Viterbi algorithm.
#' 
#' @param x Data vector the model was fitted to.
#' @param Gamma Transition probability matrix (full matrix, not only off diagonal entries).
#' @param delta Stationary distribution.
#' @param autocor default 0, Autocorrelation vector, in suitable form.
#' @param mu Optimized vector of the mu parameter in the gamma distribution.
#' @param sigma Optimized vector of the mu parameter in the gamma distribution.
#' @param N Number of states.
#' @param p Degree of autocorrelation (0 - no autocorrelation).
#' 
#' @return Estimated states using Viterbi.
#' 
#' @export
#' @rdname viterbi_arp_mu
viterbi_arp_mu <-function(x, Gamma, delta, autocor=0, 
                          mu, sigma, N, p){
  n <- length(x)
  allprobs <- matrix(1,n,N)
  
  if (p==0){
    ind <- which(!is.na(x))
    
    for (j in 1:N){
      allprobs[ind,j] <- dgamma(x[ind],
                                shape=mu[j]^2/sigma[j]^2,
                                scale=sigma[j]^2/mu[j])
    }
  } else{ # p!=0, autocorrelation
    ind <- which(!is.na(x))[-c(1:p)] # change: we omit first step 
    # in order to always have the step in t-1
    
    autocor <- matrix(autocor, ncol=p, byrow=TRUE) # autocorrelation matrix for easier handling later on
    
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


#' Global decoding for AR(p)-gamma HMM using Viterbi
#' 
#' Global decoding for an AR(p)-gamma HMM with autocorrelation in the parameters \eqn{\mu}
#' and \eqn{\sigma}
#' using the Viterbi algorithm.
#' 
#' @param x Data vector the model was fitted to.
#' @param Gamma Transition probability matrix (full matrix, not only off diagonal entries).
#' @param delta Stationary distribution.
#' @param autocor default 0, Autocorrelation vector, in suitable form.
#' @param mu Optimized vector of the mu parameter in the gamma distribution.
#' @param sigma Optimized vector of the mu parameter in the gamma distribution.
#' @param N Number of states.
#' @param p Degree of autocorrelation (0 - no autocorrelation).
#' 
#' @return Estimated states using Viterbi.
#' 
#' @export
#' @rdname viterbi_gamma_arp
viterbi_gamma_arp <-function(x, Gamma, delta, autocor=0, 
                             mu, sigma, N, p){
  n <- length(x)
  allprobs <- matrix(1,n,N)
  cv <- sigma/mu # coefficient of variance
  
  if (p==0){
    ind <- which(!is.na(x))
    
    for (j in 1:N){
      allprobs[ind,j] <- dgamma(x[ind],
                                shape=mu[j]^2/sigma[j]^2,
                                scale=sigma[j]^2/mu[j])
    }
  } else{ # p!=0, autocorrelation
    ind <- which(!is.na(x))[-c(1:p)] # change: we omit first step 
    # in order to always have the step in t-1
    
    autocor <- matrix(autocor, ncol=p, byrow=TRUE) # autocorrelation matrix for easier handling later on
    
    autocor_ind <- matrix(NA,nrow=length(ind),ncol=p) # matrix for indices of autocor data
    for (i in 1:p){
      autocor_ind[,i] <- ind-p+i-1
    }
    autocor_ind <- apply(autocor_ind, 2, function(a)x[a]) # substitute indices with values
    
    for (j in 1:N){
      mu_auto <- c(rep(NA,p), # AR(p)
                   ((1-sum(autocor[j,]))*mu[j] + 
                      as.vector(autocor_ind%*%autocor[j,]))) # matmul of values with autocor coefficient
      sigma_auto <- cv[j]*mu_auto # calculate sigma using ccv
      
      allprobs[ind,j] <- dgamma(x[ind],
                                shape=mu_auto[ind]^2/sigma_auto[ind]^2,
                                scale=sigma_auto[ind]^2/mu_auto[ind])
    }
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



#' Global decoding for AR(p)-von Mises HMM using Viterbi
#' 
#' Global decoding for an AR(p)-von Mises HMM with autocorrelation in the parameter \eqn{\mu}
#' using the Viterbi algorithm.
#' 
#' @param x Data vector the model was fitted to.
#' @param Gamma Transition probability matrix (full matrix, not only off diagonal entries).
#' @param delta Stationary distribution.
#' @param autocor default 0, Autocorrelation vector, in suitable form.
#' @param mu Optimized vector of the mu parameter in the von Mises distribution.
#' @param kappa Optimized vector of the kappa parameter in the von Mises distribution.
#' @param N Number of states.
#' @param p Degree of autocorrelation (0 - no autocorrelation).
#' 
#' @return Estimated states using Viterbi.
#' 
#' @export
#' @rdname viterbi_vonMises_arp
viterbi_vonMises_arp <-function(x, Gamma, delta, autocor=0, 
                                mu, kappa, N, p){
  n <- length(x)
  allprobs <- matrix(1,n,N)
  
  if (p==0){
    ind <- which(!is.na(x))
    
    for (j in 1:N){
      allprobs[ind,j] <- dvm(x[ind],mu[j],kappa[j])
    }
  } else{ # p!=0, autocorrelation
    ind <- which(!is.na(x))[-c(1:p)] # change: we omit first step 
    # in order to always have the step in t-1
    
    autocor <- matrix(autocor, ncol=p, byrow=TRUE) # autocorrelation matrix for easier handling later on
    
    autocor_ind <- matrix(NA,nrow=length(ind),ncol=p) # matrix for indices of autocor data
    for (i in 1:p){
      autocor_ind[,i] <- ind-p+i-1
    }
    autocor_ind <- apply(autocor_ind, 2, function(a)x[a]) # substitute indices with values
    
    for (j in 1:N){
      mu_auto <- c(rep(NA,p), # AR(p)
                   ((1-sum(autocor[j,]))*mu[j] + 
                      as.vector(autocor_ind%*%autocor[j,]))) # matmul of values with autocor coefficient
      
      allprobs[ind,j] <- dvm(x[ind],mu_auto[ind],kappa[j])
    }
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

