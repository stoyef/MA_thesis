# 2022-06-03 
# Functions that decode the states of different HMMs


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




