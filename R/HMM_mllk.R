# 2022-06-03
# Functions to calculate the negative log-Likelihood of different HMMs



###### general function for arbitrary distribution


#' Compute negative (penalized) log-likelihood of an AR(p)-HMM
#'
#' Compute the negative log-likelihood, using one or several specified distributions.
#' The distributions have to be specified by their commonly known abbreviation in R, 
#' e.g. one of ['gamma', 'vm', 'pois', 'binom',...].
#' The named list of parameters (one value for each parameter and for each state) have to be
#' in suitable form, i.e. a vector of the working parameters.
#' In the likelihood computation, contemporaneous independence is assumed.
#' Includes an optional penalization term \eqn{\lambda} for parameter selection of \eqn{p}.
#' 
#' @param theta.star Vector of parameters in the following order: 1) Off-diagonal entries of TPM,
#'                   2) Distribution parameters for each state (each distribution at a time, i.e. 
#'                   first all parameters of dist1, then all parameters of dist2 etc.), 
#'                   3) Autoregression parameters (each distribution at a time, i.e. 
#'                   first all autoregression parameters of dist1, then all autocorrelation 
#'                   parameters of dist2 etc.) -> degree parameters for each state of each variable.
#' @param dists Vector containing abbreviated names (in R-jargon) of the distributions 
#'              to be considered in the likelihood computation.
#' @param x Data vector or matrix for which the negative log-likelihood should be computed.
#' @param N Number of states.
#' @param p_auto Vector of autoregression degrees, one value for each state of each variable 
#'               (in case of penalization choose upper bound of number of parameters).
#' @param lambda Complexity penalty (≥0) for autoregression parameters \eqn{\phi}.(default: 0 no penalization).
#' @param scale_kappa Default 1, Scaling factor for kappa to avoid numerical issues in optimization for large kappa.
#' @param zero_inf Default FALSE, indicates if the gamma distributed variables should incorporate zero-inflation.
#' @param alt_data Default NULL, provide data here if variable name x is already taken in wrapper function.
#' 
#' @return Negative (penalized) log-likelihood.
#' 
#' @export
#' @rdname mllk
mllk <- function(theta.star, dists, x, N, p_auto, lambda=0, scale_kappa=1, zero_inf=FALSE, alt_data=NULL){
  
  if (!is.null(alt_data)){ # safeguard if x can't be given as argument to mllk 
                           # (problem in computation of hessian)
    x=alt_data
  }
  # First: Working to natural parameters, list structure for better handling
  # We currently only use distributions with 2 parameters, once we use Poisson distribution etc, we need a re-write
  all_params = unstarize(theta.star=theta.star, N=N, p=p_auto, dists=dists, scale_kappa = scale_kappa, zero_inf = zero_inf)
  Gamma = all_params$Gamma
  delta = all_params$delta
  autocor = all_params$autocor
  params = all_params$params
  
  # transform data to matrix, if necessary
  if (is.vector(x)) x <- matrix(x, nrow=length(x))
  
  # allprobs calculation in separate function allprobs
  allprobs = allprobs(x=x, dists=dists, autocor = autocor, params=params, N=N, p=p_auto)
  
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:dim(x)[1]){
    foo <- phi%*%Gamma%*%diag(allprobs[t,]) # here it can happen that 
    # allprobs[t,] = c(0,0) due to numeric issues. 
    # Then the function fails
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  ## Lasso penalty for autoregression parameters
  if (lambda>0){
    sum_auto = 0
    auto = c()
    for (dist in 1:length(dists)){
      for (state in 1:N){
        #sum_auto = sum_auto + sum(abs(autocor[[dist]][[state]]),na.rm = T)
        auto = c(auto, autocor[[dist]][[state]])
      }
    }
    sum_auto = sum(sqrt((auto+1e-10)^2), na.rm = T) # Approx wie in Marius Paper
    l <- l - lambda * sum_auto
  }
  return(-l)
}




#' Matrix of all probabilities
#' 
#' Calculate matrix of all probabilities for an AR(p)-HMM with or without within-state autoregression.
#' 
#' @param x Data matrix the model was fitted to.
#' @param dists Vector of distributions in R-jargon.
#' @param autocor list of autoregression vectors, respecting order of dists, one vector for each state.
#' @param params List of optimized parameters, returned by fitting the HMM.
#' @param N Number of states.
#' @param p Vector of degree of autoregression for each distribution (0 = no autoregression).
#' 
#' @return Matrix of all probabilities.
#' 
#' @export
#' @rdname allprobs
allprobs <- function(x, dists, autocor=0, params, N, p){
  if (is.vector(x)) x <- matrix(x, nrow=length(x))
  allprobs <- matrix(1,dim(x)[1],N)
  
  for (dist in 1:length(dists)){
    if (any(p[((dist-1)*N) + 1:N]>0)){
      
      for (j in 1:N){ 
        if (p[(dist-1)*N + j] > 0){ # check if current state has autoregression
          
          ind <- which(!is.na(x[,dist]))[-c(1:p[((dist-1)*N) + j])] 
          autocor_ind <- matrix(NA,nrow=length(ind),ncol=p[((dist-1)*N) + j]) 
          
          for (i in 1:p[((dist-1)*N) + j]){
            autocor_ind[,i] <- ind-p[((dist-1)*N) + j]+i-1
          }
          
          x_wo_na = x
          x_wo_na[which(is.na(x_wo_na[,dist])),dist] = mean(x_wo_na[,dist],na.rm = TRUE)
          autocor_ind <- apply(autocor_ind, 2, function(a)x_wo_na[a,dist]) 
          
          theta_j <- params[[dist]]
          for (i in names(theta_j)) theta_j[i][[1]] = theta_j[i][[1]][j] 
          
          allprobs[ind,j] <- allprobs[ind,j] * 
            match.fun(paste('dens_', dists[dist], sep=""))(x[ind,dist], theta_j, 
                                                           autocor_ind, 
                                                           autocor[[dist]][[j]],
                                                           p[((dist-1)*N) + j])
        } else{ # current state has no autoregression
          ind <- which(!is.na(x[,dist]))
          for (j in 1:N){
            
            theta_j <- params[[dist]]
            for (i in names(theta_j)) theta_j[i][[1]] = theta_j[i][[1]][j] 
            
            allprobs[ind,j] <- allprobs[ind,j] * 
              match.fun(paste('dens_', dists[dist], sep=""))(x[ind,dist], theta_j)
          }
        }
      }
    } else{ 
      
      ind <- which(!is.na(x[,dist]))
      
      for (j in 1:N){
        theta_j <- params[[dist]]
        # theta_j consists the parameters of state j for each parameter in theta$params
        for (i in names(theta_j)) theta_j[i][[1]] = theta_j[i][[1]][j] 
        
        allprobs[ind,j] <- allprobs[ind,j] * 
          match.fun(paste('dens_', dists[dist], sep=""))(x[ind,dist], theta_j)
      }
    }
  }
  return(allprobs)
}



