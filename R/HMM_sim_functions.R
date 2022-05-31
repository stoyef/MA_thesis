# 2022-05-30
# Functions for simulating and evaluating HMMs with or without autocorrelation
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # RStudio way


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


#' Simulate data from a Gamma HMM with AR(p) structure
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
#' @rdname sample_hmm_arp
sample_hmm_arp <- function(n_samples, delta, Gamma, N, mu, sigma, autocor, p){
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
#' 
#' @return Minus Log-Likelihood
#' 
#' @export
#' @rdname mllk_hmm
mllk_hmm <-function(theta.star,x,N){
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

#' Calculate negative Log-Likelihood of gamma HMM
#' 
#' Calculate the negative Log-Likelihood of an \emph{AR(p)}-gamma HMM 
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
#' @rdname mllk_arp
mllk_arp <-function(theta.star,x,N,p){
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  
  autocor <- plogis(theta.star[(N-1)*N+1:(p*N)])
  autocor <- matrix(autocor, ncol=p, byrow=TRUE) # matrix for easier handling later on
  mu <- exp(theta.star[(N-1)*N+(p*N+1):(p*N+N)])
  sigma <- exp(theta.star[(N-1)*N+p*N+N+1:N])
  
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


#' Fit an AR(p)-gamma HMM to data
#' 
#' Fit an AR(p)-gamma HMM to data using a specified function to compute the negative 
#' Log-Likelihood. This function gets minimized using the function \code{optim}.
#' It returns the estimated parameters of the fitted model.
#' 
#' @param mllk Negative log-Likelihood function that should be minimized.
#' @param data Data that should be fitted using the HMM.
#' @param theta.star Unconstrained parameter vector (has to be provided in suitable form).
#' @param N Number of states.
#' @param p Degree of autocorrelation, 0 is possible for no autocorrelation.
#' 
#' @return List, containing minimal value of negative log-Likelihood, Gamma, 
#'         delta, (autocorrelation, depending on degree), mu, sigma.
#' 
#' @export
#' @rdname fit_arp_model
fit_arp_model <- function(mllk, data, theta.star, N, p){
  if (p>0){
    mod <- optim(par=theta.star, fn=mllk, method='L-BFGS-B',
                 N=N,p=p,x=data)
  } else{
    mod <- optim(par=theta.star, fn=mllk, method='L-BFGS-B',
                 N=N,x=data)
  }
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(mod$par[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  if (p>0){
    # autocorrelation
    autocor <- plogis(mod$par[(N-1)*N+1:(p*N)])
    # mu and sigma 
    mu <- exp(mod$par[(N-1)*N+(p*N+1):(p*N+N)])
    sigma <- exp(mod$par[(N-1)*N+p*N+N+1:N])
    ret <- list(mod$value, Gamma, delta, autocor, mu, sigma)
    names(ret) <- c('mllk', 'Gamma', 'delta', 'autocorrelation', 'mu', 'sigma')
    return(ret)
  } else{
    # mu and sigma 
    mu <- exp(mod$par[(N-1)*N+1:N])
    sigma <- exp(mod$par[(N-1)*N+N+1:N])
    ret <- list(mod$value, Gamma, delta, mu, sigma)
    names(ret) <- c('mllk', 'Gamma', 'delta', 'mu', 'sigma')
    return(ret)
  }
}


#' Calculate AIC for fitted gamma HMM
#'
#' Calculate the AIC criterion for a fitted gamma HMM 
#' (with or without autocorrelation).
#' 
#' @param mllk minimum negative Log-Likelihood, output of \code{fit_gamma_hmm}.
#' @param N Number of states of the HMM.
#' @param p Degree of autocorrelation of the HMM (0 - no autocorrelation).
#' 
#' @return AIC of fitted model.
#' 
#' @export
#' @rdname AIC_gamma_HMM
AIC_gamma_HMM <- function(mllk, N, p){
  n_params <- N*(N-1) + p*N + 2*N # TPM+autocor+mu+sigma
  ret <- 2*mllk + 2*n_params # AIC (logLike is already negative)
  return(ret)
}


#' Calculate BIC for fitted gamma HMM
#'
#' Calculate the BIC criterion for a fitted gamma HMM 
#' (with or without autocorrelation).
#' 
#' @param mllk minimum negative Log-Likelihood, output of \code{fit_gamma_hmm}.
#' @param N Number of states of the HMM.
#' @param p Degree of autocorrelation of the HMM (0 - no autocorrelation).
#' @param data Data vector the data was fitted to, to determine \eqn{\log(n)}.
#' 
#' @return BIC of fitted model.
#' 
#' @export
#' @rdname BIC_gamma_HMM
BIC_gamma_HMM <- function(mllk, N, p, data){
  n_params <- N*(N-1) + p*N + 2*N # TPM+autocor+mu+sigma
  ret <- 2*mllk + log(length(data))*n_params # BIC (logLike is already negative)
  return(ret)
}


#' Global decoding for AR(p)-gamma HMM using Viterbi
#' 
#' Global decoding for an AR(p)-gamma HMM using the Viterbi algorithm.
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
#' @rdname viterbi_arp
viterbi_arp <-function(x, Gamma, delta, autocor=0, 
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


#' Plot states decoded by Viterbi algorithm
#' 
#' This function plots the states decoded by the Viterbi algorithm.
#' 
#' @param states Vector of states (natural numbers, range 1-N).
#' @param names optional, specifies the names of the states. If "none", state 1-N will be used.
#' @param title bool, indicator if the plot should habe a title.
#' 
#' @export
#' @rdname plot_states
plot_states <- function(states, names="none", title=TRUE){
  N <- length(unique(states))
  require(RColorBrewer)
  pal <- brewer.pal(N,'Dark2')
  
  if (title){
  plot(states, pch=19, bty='n', main='Decoded states using Viterbi',
       xlab='Index', ylab='State', col=pal[states])
  } else{
    plot(states, pch=19, bty='n', main='',
         xlab='Index', ylab='State', col=pal[states])
    }
  if (names[1]=="none"){
    legend(length(states)/1.5, 1.8, paste("State",1:N), col=pal,lwd=1.5,bty='n')
  } else{
    legend(length(states)/1.5, 1.8, names, col=pal,lwd=1.5,bty='n')
  }
}


#' Density plot of estimated gamma distributions
#'
#' Plots the estimted gamma distributions along with a histogram of the values.
#' 
#' @param data Input data for histogram
#' @param mu Input vector for parameter mu
#' @param sigma Input vector for parameter sigma
#' @param delta Input vector for parameter delta
#' @param title optional, title of the plot, if "none", then no title
#' 
#' @export
#' @rdname plot_fitted_gamma_dist
plot_fitted_gamma_dist <- function(data, mu, sigma, delta, title="none"){
  N <- length(mu)
  require(RColorBrewer)
  pal <- brewer.pal(N+1, 'Dark2')
  gams <- matrix(NA, nrow=10000,ncol=N)
  x <- seq(min(data),max(data),length=10000)
  for (i in 1:N){
    gams[,i] = delta[i]*dgamma(x, shape=mu[i]^2/sigma[i]^2, scale=sigma[i]^2/mu[i])
  }
  total_dist <- apply(gams,1,sum)
  
  if (title=="none"){
    hist(data, breaks=length(data)/25, probability = TRUE,
         main="", xlab="x", ylim=c(0,1.1*max(total_dist)))
  }else{
  hist(data, breaks=length(data)/25, probability = TRUE,
       main=title, xlab="x", ylim=c(0,1.1*max(total_dist)))
  }
  for (i in 1:N){
    lines(x,gams[,i],col=pal[i], lwd=2)
  }
  lines(x,total_dist,col=pal[N+1],lwd=2)
  legend('topright', c(paste("State",1:N),"Total"), bty='n', lwd=2,
         col=pal)
}


#' Plot data 
#'
#' Plot a data vector
#' 
#' @param data Data vector to be plotted
#' @param name optional, y-axis label, if "none" than "data"
#' @param title optional, title of the plot, if "none", then no title
#' 
#' @export
#' @rdname plot_data
plot_data <- function(data, name="none", title="none"){
  if (name=="none"){
    plot(data, typ='l',xlab="Index", ylab="Data",bty='n')
  } else{
    plot(data, typ='l',xlab="Index", ylab=name,bty='n')
  }
  if (title=="none"){
    title(main="")
  } else{
    title(main=title)
  }
  
}


#' Simulate gamma HMMs with and without autocorrelation
#'
#' Top to bottom wrapper function: Simulate from a specified HMM, fit a specified
#' HMM to the resulting data and return the output. Optionally: Also return the 
#' decoded states using the Viterbi algorithm. Optionally: Also plot the resulting 
#' density functions of the weighted gamma distributions.
#' 
#' @param model_sim Int, Form of simulated data: Degree of autocorrelation 
#'                       (0 - no autocorrelation).
#' @param model_fit Int, Form of fitted model: Degree of autocorrelation 
#'                       (0 - no autocorrelation).
#' @param N_sim Number of states of the HMM used for data generation.
#' @param N_fit Number of states of the HMM used for data generation.
#' @param n_samples Number of samples generated.
#' @param Gamma_sim Full transition probability matrix of simulated data.
#' @param delta_sim Initial distribution of simulated data.
#' @param mu_sim Parameter vector for mu of simulated data.
#' @param sigma_sim Parameter vector for sigma of simulated data.
#' @param autocor_sim Vector of autocorrelation coefficients of simulated data (if there are any).
#' @param estimate_states Bool, determines if states are estimated and returned
#'                        using Viterbi.
#' @param plot_it Bool, determines if resulting densities are plotted.    
#' 
#' @return List of Fitted model and its parameters (and optional the decoded states).
#' 
#' @export
#' @rdname gamma_simulation            
#' 
gamma_simulation <- function(model_sim, model_fit, N_sim, N_fit, n_samples, 
                             Gamma_sim, delta_sim, mu_sim, sigma_sim, autocor_sim=0,
                             estimate_states=TRUE,plot_it=TRUE){
  
  if (model_sim>0){
    ## prepare autocor for sampling -> matrix
    autocor_matrix <- matrix(autocor_sim, ncol=model_sim, byrow=TRUE) # matrix for easier handling later on
  }
  
  ## simulated data -> contained in simulated_data$data
  if (model_sim==0){
    simulated_data <- sample_hmm_normal(n_samples, delta_sim, Gamma_sim, N_sim,
                                        mu_sim, sigma_sim)
  } else if (model_sim==1){
    simulated_data <- sample_hmm_ar1(n_samples, delta_sim, Gamma_sim, N_sim,
                                           mu_sim, sigma_sim, autocor_matrix)
  } else if (model_sim>1){
    simulated_data <- sample_hmm_arp(n_samples, delta_sim, Gamma_sim, N_sim,
                                     mu_sim, sigma_sim, autocor_matrix,p=model_sim)
  } else{
    return("Wrong input for parameter model_sim.")
  }
  
  # Parameters for model fitting
  theta = c(
    rep(-2,N_fit*(N_fit-1)), # TPM
    rep(0.1,N_fit*model_fit)/(N_fit*model_fit+1), # autocor, n_states*p, regularization to avoid sum > 1
    quantile(simulated_data$data, seq(0,1,length=N_fit)), # mu 
    rep(sd(simulated_data$data),N_fit) # sigma
  )
  
  if (model_fit==0){
    theta.star = c(
      theta[1:(N_fit*(N_fit-1))],
      log(theta[N_fit*(N_fit-1)+1:(2*N_fit)])
    )
  } else{
    theta.star = c(
      theta[1:(N_fit*(N_fit-1))],
      qlogis(theta[(N_fit*(N_fit-1))+1:(N_fit*model_fit)]),
      log(theta[N_fit*(N_fit-1)+(N_fit*model_fit)+1:(2*N_fit)])
    )
  }

  if (model_fit==0){ # no autocorrelation in fitted model
    fitted_model <- fit_arp_model(mllk_hmm, simulated_data$data,
                                  theta.star, N=N_fit, p=model_fit)
  } else if (model_fit>0){ # there is autocorrelation in fitted model
    fitted_model <- fit_arp_model(mllk_arp, simulated_data$data,
                                  theta.star, N=N_fit, p=model_fit)
  } else{
    return("Wrong input for parameter model_fit.")
  }
  
  # Plot, if wanted
  if (plot_it){
    plot_fitted_gamma_dist(simulated_data$data, fitted_model$mu, fitted_model$sigma,
                           fitted_model$delta)
  }
  
  # Viterbi, if wanted
  if (estimate_states){
    estimated_states <- viterbi_arp(simulated_data$data, fitted_model$Gamma,
                                    fitted_model$delta, fitted_model$autocor,
                                    fitted_model$mu, fitted_model$sigma,
                                    N_fit, model_fit)
    
    ret <- list(simulated_data, fitted_model, estimated_states)
    names(ret) <- c('simulated_model','fitted_model', 'viterbi_states')
    return(ret)
  } else{
    ret <- list(simulated_data, fitted_model)
    names(ret) <- c('simulated_model','fitted_model')
    return(ret)
  }
}

