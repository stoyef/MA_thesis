# 2022-06-03 
# Functions that simulate data from different HMMs




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











################################################################################
################################################################################
## Below are only old functions that are not in use currently
## They probably can be deleted in the future.



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



#' Simulate data from an HMM with AR(p) structure
#' 
#' Simulate data from an HMM with AR(p) structure. 
#' Different distributions can be specified in \code{dists} (uni- and multivariate).
#' 
#'
#' 
#' @param n_samples Number of samples to generate.
#' @param delta Initial distribution of the Markov chain.
#' @param Gamma Transition probability matrix of the Markov chain.
#' @param N Number of states.
#' @param params Parameter vector for the different distributions. 
#'.              Has to respect the order specified in \code{dists}.
#' @param autocor List of parameter matrix (Dimensions for each matrix: \eqn{n\times p}) for the autocorrelation coefficients. 
#'             Has to match p, in the order \eqn{\phi_{t-p},\dots,\phi_{t-1}}
#'             where \eqn{\phi} is the vector of autocorrelation coefficients
#'             for one specific time lag (one value for each state).
#'             0, if no autocorrelation. Has to respect the order specified in \code{dists}.
#' @param p Vector of degree of autocorrelation for each distribution, 0=no autocorrelation.
#' @param dists Vector containing abbreviated names (in R-jargon) of the distributions 
#'              to be considered in the Likelihood computation.
#'               
#' @return List of states and data of the HMM.
#' 
#' @export
#' @rdname sample_arp
sample_arp <- function(n_samples, delta, Gamma, N, params, autocor, p, dists){
  
  states <- rep(NA,n_samples)
  data <- matrix(NA,nrow=n_samples,ncol=length(dists))
  # calculate cv in case we need them
  cv = rep(NA, length(dists)*N)
  s = N
  m = 0
  for (i in 1:length(dists)){
    cv[(i-1)*N+1:N] = params[s+1:N] / params[m+1:N]
    s = s + 2*N
    m = m + 2*N 
  }
  
  # first p data points given, no autocorrelation -> just like normal HMM
  # we have to be careful with rvm, it uses 0-2*pi instead of -pi-+pi
  # currently only works for distributions with 2 parameters!
  states[1] <- sample(1:N, 1, prob = delta)
  for (dist in 1:length(dists)){
    param1=params[(N*(dist-1)*N)+1:N] # first parameter
    param2=params[(N*(dist-1)*N)+N+1:N] # second parameter
    # choose distribution and sample
    data[1,dist] <- match.fun(paste('sample_', dists[dist], sep=""))(1, param1[states[1]], param2[states[1]])
  }
  if (p>1){
    for (t in 2:p){
      states[t] <- sample(1:N, 1, prob = Gamma[states[t-1],])
      for (dist in 1:length(dists)){
        param1=params[(N*(dist-1)*N)+1:N] # first parameter
        param2=params[(N*(dist-1)*N)+N+1:N] # second parameter
        # choose distribution and sample
        data[t,dist] <- match.fun(paste('sample_', dists[dist], sep=""))(1, param1[states[t]], param2[states[t]])
      }
    }
  }
  if (p==0){ # no autocorrelation
    for (t in 2:n_samples){
      states[t] <- sample(1:N, 1, prob = Gamma[states[t-1],])
      for (dist in 1:length(dists)){
        param1=params[(N*(dist-1)*N)+1:N] # first parameter
        param2=params[(N*(dist-1)*N)+N+1:N] # second parameter
        # choose distribution and sample
        data[t,dist] <- match.fun(paste('sample_', dists[dist], sep=""))(1, param1[states[t]], param2[states[t]])
      }
    }
  } else{ # autocorrelation
    for (t in (p+1):n_samples){
      states[t] <- sample(1:N, 1, prob=Gamma[states[t-1],])
      for (dist in 1:length(dists)){
        param1=params[(N*(dist-1)*N)+1:N] # first parameter
        param2=params[(N*(dist-1)*N)+N+1:N] # second parameter
        ar_matrix <- autocor[[dist]]
        param1_ar <- sum(ar_matrix[states[t],]*data[(t-p):(t-1),dist]) + 
          (1-sum(ar_matrix[states[t],]))*param1[states[t]]
        
        if (dists[dist]=='gamma'){ # respect ccv
          param2_ar=cv[N*(dist-1)+1:N] * param1_ar
          param2_ar = param2_ar[states[t]]
        } else{
          param2_ar=param2[states[t]]
        }
        
        # choose distribution and sample
        data[t,dist] <- match.fun(paste('sample_', dists[dist], sep=""))(1, param1_ar, param2_ar)
      }
    }
  }
  
  ret <- list(states, data)
  names(ret) <- c('states','data')
  return(ret)
}






