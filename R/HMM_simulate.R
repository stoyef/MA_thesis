# 2022-06-03 
# Functions that simulate data from different HMMs




#' Simulate data from an AR(p)-HMM
#' 
#' Simulate data from an AR(p)-HMM with AR(p).
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
#' @param autocor List of parameter matrix (Dimensions for each matrix: \eqn{n\times p}) for the autoregression parameters. 
#'             Has to match p, in the order \eqn{\phi_{t-p},\dots,\phi_{t-1}}
#'             where \eqn{\phi} is the vector of autoregression parameters
#'             for one specific time lag (one value for each state).
#'             0, if no autoregression Has to respect the order specified in \code{dists}.
#' @param p Vector of degree of autoregression for each distribution, 0 = no autoregression
#' @param dists Vector containing abbreviated names (in R-jargon) of the distributions 
#'              to be considered in the likelihood computation.
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
    param1=params[(2*(dist-1)*N)+1:N] # first parameter
    param2=params[(2*(dist-1)*N)+N+1:N] # second parameter
    # choose distribution and sample
    data[1,dist] <- match.fun(paste('sample_', dists[dist], sep=""))(1, param1[states[1]], param2[states[1]])
  }
  if (any(p>1)){
    for (t in 2:max(p)){ # we assume the data as given for all distributions until time point p
      states[t] <- sample(1:N, 1, prob = Gamma[states[t-1],])
      for (dist in 1:length(dists)){
        param1=params[(2*(dist-1)*N)+1:N] # first parameter
        param2=params[(2*(dist-1)*N)+N+1:N] # second parameter
        # choose distribution and sample
        data[t,dist] <- match.fun(paste('sample_', dists[dist], sep=""))(1, param1[states[t]], param2[states[t]])
      }
    }
  }
  if (sum(p==0) == length(p)){ # no autocorrelation
    for (t in 2:n_samples){
      states[t] <- sample(1:N, 1, prob = Gamma[states[t-1],])
      for (dist in 1:length(dists)){
        param1=params[(2*(dist-1)*N)+1:N] # first parameter
        param2=params[(2*(dist-1)*N)+N+1:N] # second parameter
        # choose distribution and sample
        data[t,dist] <- match.fun(paste('sample_', dists[dist], sep=""))(1, param1[states[t]], param2[states[t]])
      }
    }
  } else{ # autocorrelation
    for (t in (max(p)+1):n_samples){
      states[t] <- sample(1:N, 1, prob=Gamma[states[t-1],])
      for (dist in 1:length(dists)){
        param1=params[(2*(dist-1)*N)+1:N] # first parameter
        param2=params[(2*(dist-1)*N)+N+1:N] # second parameter
        ar_matrix <- autocor[[dist]]
        
        if (dists[dist]=='vm'){ # compute mean value on a circle
          param1_ar_auto <- sum(ar_matrix[states[t],]*exp(complex(imaginary=data[(t-p[dist]):(t-1),dist])))
          param1_ar_fixed <- (1-sum(ar_matrix[states[t],]))*exp(complex(imaginary=param1[states[t]]))
          param1_ar <- Arg(param1_ar_auto + param1_ar_fixed)
        } else{
          param1_ar <- sum(ar_matrix[states[t],]*data[(t-p[dist]):(t-1),dist]) + 
            (1-sum(ar_matrix[states[t],]))*param1[states[t]]
        }
        
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






