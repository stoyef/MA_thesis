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



#' Global decoding for AR(p) HMM using Viterbi
#' 
#' Global decoding for an AR(p) HMM with or without autocorrelation 
#' using the Viterbi algorithm.
#' 
#' @param x Data matrix the model was fitted to.
#' @param Gamma Transition probability matrix (full matrix, not only off diagonal entries).
#' @param delta Stationary distribution.
#' @param dists Vector of distributions in R-jargon.
#' @param autocor list of autocorrelation vectors, respecting order of dists.
#' @param params List of optimized parameters, returned by fitting the HMM.
#' @param N Number of states.
#' @param p Vector of degree of autocorrelation for each distribution (0 = no autocorrelation).
#' 
#' @return Estimated states using Viterbi.
#' 
#' @export
#' @rdname viterbi_arp
viterbi_arp <-function(x, Gamma, delta, dists, autocor=0, 
                             params, N, p){
  n <- dim(x)[1]
  allprobs <- matrix(1,n,N)
  
  for (dist in 1:length(dists)){
  if (p[dist]==0){
      ind <- which(!is.na(x[,dist]))
      
      for (j in 1:N){
        params_j <- params[[dist]]
        # params_j consists the parameters of state j for each parameter in params[[dist]]
        for (i in names(params_j)) params_j[i][[1]] = params_j[i][[1]][j] 
        allprobs[ind,j] <- allprobs[ind,j] * match.fun(paste('dens_', dists[dist], sep=""))(x[ind,dist], 
                                                                                            params_j)
      }
    } else{ # p!=0, autocorrelation
      
      ind <- which(!is.na(x[,dist]))[-c(1:p[dist])] # change: we omit first step 
      # in order to always have the step in t-1
      
      autocor_m <- matrix(autocor[[dist]], ncol=p[dist], byrow=TRUE) # autocorrelation matrix for easier handling later on
      
      autocor_ind <- matrix(NA,nrow=length(ind),ncol=p[dist]) # matrix for indices of autocor data
      for (i in 1:p[dist]){
        autocor_ind[,i] <- ind-p[dist]+i-1
      }
      autocor_ind <- apply(autocor_ind, 2, function(a)x[a,dist]) # substitute indices with values
      
      for (j in 1:N){
        params_j <- params[[dist]]
        # params_j consists the parameters of state j for each parameter in params[[dist]]
        for (i in names(params_j)) params_j[i][[1]] = params_j[i][[1]][j] 
        
        allprobs[ind,j] <- allprobs[ind,j] * 
          match.fun(paste('dens_', dists[dist], sep=""))(x[ind,dist], 
                                                         params_j,
                                                         autocor_ind=autocor_ind,
                                                         autocor=autocor_m[j,],
                                                         p=p[dist])
      }
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



