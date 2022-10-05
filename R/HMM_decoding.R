# 2022-06-03 
# Functions that decode the states of different HMMs


#' Matrix of all probabilities
#' 
#' Calculate matrix of all probabilities for an AR(p)-HMM with or without within-state autoregression.
#' 
#' @param x Data matrix the model was fitted to.
#' @param dists Vector of distributions in R-jargon.
#' @param autocor list of autoregression vectors, respecting order of dists.
#' @param params List of optimized parameters, returned by fitting the HMM.
#' @param N Number of states.
#' @param p Vector of degree of autoregression for each distribution (0 = no autoregression).
#' 
#' @return Matrix of all probabilities.
#' 
#' @export
#' @rdname allprobs
allprobs <- function(x, dists, autocor=0, params, N, p){
  if (length(dists)>1){
    n <- dim(x)[1]
  } else{
    n <- length(x)
  }
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
    } else{ # p!=0, autoregression
      
      if (length(dists)>1){
        ind <- which(!is.na(x[,dist]))[-c(1:p[dist])] # change: we omit first step 
        # in order to always have the step in t-1
      } else{
        ind <- which(!is.na(x))[-c(1:p[dist])]
      }
      
      autocor_m <- matrix(autocor[[dist]], ncol=p[dist], byrow=TRUE) # autoregression matrix for easier handling later on
      
      autocor_ind <- matrix(NA,nrow=length(ind),ncol=p[dist]) # matrix for indices of autocor data
      for (i in 1:p[dist]){
        autocor_ind[,i] <- ind-p[dist]+i-1
      }
      
      x_wo_na = x
      x_wo_na[which(is.na(x_wo_na[,dist])),dist] = mean(x_wo_na[,dist],na.rm = TRUE)
      # replace NA values in x that are put into autocor_ind with mean value 
      # (so that the data that has NA in previous time steps does not have to be deleted)
      if (length(dists)>1){
        autocor_ind <- apply(autocor_ind, 2, function(a)x_wo_na[a,dist]) # substitute indices with values
      } else{
        autocor_ind <- apply(autocor_ind, 2, function(a)x_wo_na[a])
      }
      
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
  return(allprobs)
}


#' Global decoding for AR(p)-HMM using Viterbi
#' 
#' Global decoding for an AR(p)-HMM with or without within-state autoregression 
#' using the Viterbi algorithm.
#' 
#' @param x Data matrix the model was fitted to.
#' @param Gamma Transition probability matrix (full matrix, not only off diagonal entries).
#' @param delta Stationary distribution.
#' @param dists Vector of distributions in R-jargon.
#' @param autocor list of autoregression vectors, respecting order of dists.
#' @param params List of optimized parameters, returned by fitting the HMM.
#' @param N Number of states.
#' @param p Vector of degree of autoregression for each distribution (0 = no autocorrelation).
#' 
#' @return Estimated states using Viterbi.
#' 
#' @export
#' @rdname viterbi_arp
viterbi_arp <-function(x, Gamma, delta, dists, autocor=0, 
                             params, N, p){

  allprobs = allprobs(x=x, dists=dists, autocor=autocor, params=params, N=N, p=p)  
  if (length(dists)>1){
    n <- dim(x)[1]
  } else{
    n <- length(x)
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



