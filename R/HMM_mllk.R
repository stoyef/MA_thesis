# 2022-06-03
# Functions to calculate the negative log-Likelihood of different HMMs



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


###### general function for arbitrary distribution


#' Compute negative Log-Likelihood
#'
#' Compute the negative Log-Likelihood, using one or several specified distributions.
#' The distributions have to be specified by their commonly known abbreviation in R, 
#' e.g. one of ['gamma', 'vm', 'pois', 'binom',...].
#' The named list of parameters (one value for each parameter and for each state) have to be
#' in suitable form, i.e. they have to respect the natural parameter boundaries of the distributions.
#' In the likelihood computation, contemporaneus independence is assumed.
#' 
#' @param theta.star Named list of parameters, containing: 1) Full TPM, 2) List of 
#'              autocorrelation parameters, one vector entry for each distribution, 
#'              3) List with sublist for each distribution, containing named 
#'              parameters (N entries for each parameter of each distribution)
#' @param dists Vector containing abbreviated names (in R-jargon) of the distributions 
#'              to be considered in the Likelihood computation.
#' @param x Data vector or matrix for which the negative Log-Likelihood should be computed.
#' @param N Number of states.
#' @param p Vector of degree of autocorrelation for each distribution, 0=no autocorrelation.
#' 
#' @return Negative Log-Likelihood.
#' 
#' @export
#' @rdname mllk
mllk <- function(theta.star, dists, x, N, p){
  
  # First: Working to natural parameters, list structure for better handling
  # We currently only use distributions with 2 parameters, once we use Poisson distribution etc, we need a re-write
  
  ### TPM
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma) 
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  
  ### distribution parameters
  n_dists <- length(dists)
  n_dist_params <- 2*n_dists*N
  params = list()
  counter = N*(N-1)
  for (dist in 1:n_dists){
    if (dists[dist]=='gamma'){
      params[[dist]] = list(mu=exp(theta.star[counter+1:N]),
                          sigma=exp(theta.star[counter+N+1:N]))
    } else if (dists[dist]=='vm'){
      params[[dist]] = list(mu = Arg(theta.star[counter+1:N]+1i*theta.star[counter+N+1:N]),
                            kappa = sqrt(theta.star[counter+1:N]^2+theta.star[counter+N+1:N]^2))
      
    } else if (dists[dist]=='norm'){
      params[[dist]] = list(mu=theta.star[counter+1:N],
                            sigma=exp(theta.star[counter+N+1:N]))
    } else{
      return(paste("ERROR: The distribution", dists[dist], "is not implemented."))
    }
    counter = counter+2*N
  }
  
  ### Autocorrelation parameters
  a_params = theta.star[-c(1:counter)]
  autocor = list()
  counter = 0
  for (dist in 1:n_dists){
    autocor[[dist]] = plogis(a_params[counter+1:(p[dist]*N)])
    counter = counter+p[dist]*N
  }
  
  
  # transform data to matrix, if necessary
  if (is.vector(x)) x <- matrix(x, nrow=length(x))
  allprobs <- matrix(1,dim(x)[1],N)
  
  for (dist in 1:length(dists)){ # for each distribution (= column of x) to consider

    if (p[dist]>0){
      autocor_m <- matrix(autocor[[dist]], ncol=p[dist], byrow=TRUE) # matrix for easier handling later on
      
      ind <- which(!is.na(x[,dist]))[-c(1:p[dist])] # change: we omit first p steps 
      # in order to always have the step in t-p
      autocor_ind <- matrix(NA,nrow=length(ind),ncol=p[dist]) # matrix for indices of autocor data
      
      for (i in 1:p[dist]){
        autocor_ind[,i] <- ind-p[dist]+i-1
      }
      autocor_ind <- apply(autocor_ind, 2, function(a)x[a,dist]) # substitute indices with values
      
      for (j in 1:N){
        # here comes the autocorrelation! -> computed inside the dens_<...> functions, considering the 
        # autocorrelation!
        theta_j <- params[[dist]]
        # theta_j consists the parameters of state j for each parameter in theta$params
        for (i in names(theta_j)) theta_j[i][[1]] = theta_j[i][[1]][j] 
        
        # computation of the different densities is outsourced to the respective 
        # functions with name "d<name_of_distribution>, e.g. dgamma for gamma distribution
        allprobs[ind,j] <- allprobs[ind,j] * 
                              match.fun(paste('dens_', dists[dist], sep=""))(x[ind,dist], theta_j, autocor_ind, autocor_m[j,], p[dist])
        # here we have to choose mu_auto[ind], because
        # we have an individual mu for each data point
      }
    } else{
      ind <- which(!is.na(x[,dist]))
      
      for (j in 1:N){
        
        theta_j <- params[[dist]]
        # theta_j consists the parameters of state j for each parameter in theta$params
        for (i in names(theta_j)) theta_j[i][[1]] = theta_j[i][[1]][j] 
        
        allprobs[ind,j] <- allprobs[ind,j] * 
                              match.fun(paste('dens_', dists[dist], sep=""))(x[ind,dist], theta_j, autocor_ind, autocor[[dist]], p[dist])
      }
    }
      
  }
  
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:dim(x)[1]){
    foo <- phi%*%Gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}





