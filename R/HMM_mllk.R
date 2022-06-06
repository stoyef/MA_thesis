# 2022-06-03
# Functions to calculate the negative logg-Likelihood of different HMMs



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
#' @param p Degree of autocorrelation, 0=no autocorrelation
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





