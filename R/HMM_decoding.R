# 2022-06-03 
# Functions that decode the states of different HMMs



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



