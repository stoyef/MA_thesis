# 2022-06-03 
# Functions that fit different HMMs to data


#' Fit an AR(p)-HMM to data
#' 
#' Fit an AR(p)-HMM with within-state autoregression to data 
#' using a specified function to compute the negative 
#' log-likelihood. This function gets by default minimized using the function \code{optim}.
#' It returns the estimated parameters of the fitted model. The distributional assumptions for the HMM
#' can be specified using the parameter \code{dist}.
#' 
#' @param mllk Negative log-likelihood function that should be minimized.
#' @param data Data that should be fitted using the HMM.
#' @param theta.star Unconstrained parameter vector (has to be provided in suitable form).
#'                   Attention: -Inf values (resulting e.g. from supplying autocorrelation = 0) are not possible
#'                   for the optimization function. 
#' @param N Number of states.
#' @param p_auto Vector of degree of autoregression for each distribution, 0 = no autoregression, one value
#'               for every state of every variable.
#' @param dists Vector containing abbreviated names (in R-jargon) of the distributions 
#'              to be considered in the likelihood computation.
#' @param opt_fun string - Function that should be used for optimization (default: optim, one of ['optim', 'nlm', 'genoud']).
#' @param scale_kappa Default 1, Scaling factor for kappa to avoid numerical issues in optimization for large kappa in the von Mises distribution.
#' @param zero_inf Default FALSE, indicates if the gamma distributed variables should incorporate zero-inflation.
#'            
#' @return List, containing minimal value of negative log-likelihood, Gamma, 
#'         delta, (autoregression, depending on degree), mu, sigma.
#' 
#' @export
#' @rdname fit_arp_model
fit_arp_model <- function(mllk, data, theta.star, N, p_auto, dists, opt_fun='optim', scale_kappa=1, zero_inf=FALSE){
  skip = FALSE
  
  if (opt_fun=='optim'){
    tryCatch(
      # Error handling: sometimes the optim function does not work, because of some
      # error. We want to be notified that there is an error, but the execution should 
      # not be interrupted (important for for loops that use this function)
      mod <- optim(par=theta.star, fn=mllk, method='L-BFGS-B',
                   N=N,p_auto=p_auto,x=data, dists=dists, scale_kappa=scale_kappa, zero_inf=zero_inf),
      error=function(e){cat("ERROR: optim() failed. Did you supply autocorrelation parameters = 0?\n
                            Continue to next iteration.\n")
        skip <<- TRUE
      }
    )
  } else if (opt_fun=='nlm'){
    tryCatch(
      mod <- nlm(f=mllk, p=theta.star, iterlim=500,
                 N=N,p_auto=p_auto,x=data, dists=dists, scale_kappa=scale_kappa, zero_inf=zero_inf),
      error=function(e){cat("ERROR: nlm() failed.\n
                            Continue to next iteration.\n")
        skip <<- TRUE
      }
    )
  } else if (opt_fun=="ga"){
    require("rgenoud")
    tryCatch(
      mod <- genoud(fn=mllk, nvars=length(theta.star), pop.size=1000, max.generations=100,
                    starting.values = theta.star,
                    N=N,p_auto=p_auto,scale_kappa=scale_kappa,zero_inf=zero_inf,
                    dists=dists,x=data#,
                    #Domains=matrix(rep(c(-Inf,Inf),length(theta.star)),ncol=2,byrow=T)
      ),
      error=function(e){cat("ERROR: rgenoud() failed.\n
                            Continue to next iteration.\n")
        skip <<- TRUE
      }
    )
  }
  
  # leave function, if optim() didn't work
  if (skip){
    return(NA)
  }
  
  # working parameters -> natural parameters
  ## TPM
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(mod$par[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  
  # create parameter list
  params = list()
  
  counter = N*(N-1)
  # fill parameters in list structure
  for (dist in 1:length(dists)){
    if (dists[dist]=='gamma'){
      # mu and sigma 
      mu <- exp(mod$par[counter+1:N])
      sigma <- exp(mod$par[counter+N+1:N])
      params[[dist]] = list(mu=mu, sigma=sigma)
      counter = counter+2*N
    } else if (dists[dist]=='vm'){
      # mu and kappa 
      mu <- Arg((mod$par[counter+1:N])*scale_kappa+1i*(mod$par[counter+N+1:N]*scale_kappa))
      kappa <- sqrt((mod$par[counter+1:N]*scale_kappa)^2+(mod$par[counter+N+1:N]*scale_kappa)^2)
      params[[dist]] = list(mu=mu, kappa=kappa)
      counter = counter+2*N
    } else if (dists[dist]=='norm'){
      # mu and sigma 
      mu <- mod$par[counter+1:N]
      sigma <- exp(mod$par[counter+N+1:N])
      params[[dist]] = list(mu=mu, sigma=sigma)
      counter = counter+2*N
    }
  }
  
  counter_auto = 0
  # autocorrelation
  if (any(p_auto>0)){
    ac <- plogis(mod$par[counter+1:sum(p_auto)])
    autocor = list()
    for (dist in 1:length(dists)){
      autocor[[dist]] = list()
      for (state in 1:N){
        if (p_auto[(dist-1)*N+state]>0){
          autocor[[dist]][[state]] = ac[counter_auto+1:(p_auto[(dist-1)*N+state])]
          counter_auto = counter_auto+p_auto[(dist-1)*N+state]
        }
      }
    }
  }
  
  if (zero_inf){
    gamma_dists = which(dists=='gamma')
    zero_inf <- plogis(mod$par[-c(1:(counter+counter_auto))])
    counter_zero_inf = 0
    for (dist in gamma_dists){
      params[[dist]]['zero_inf'] = list(zero_inf[counter_zero_inf+1:N])
      counter_zero_inf = counter_zero_inf+N
    }
  }
  
  # AIC & BIC (mod$value is already negative log-likelihood)
  aic = 2*mod$value + 2*length(theta.star)
  bic = 2*mod$value + log(dim(data)[1])*length(theta.star)
  
  # create return object
  if (any(p_auto>0)){
    ret <- list(mod$value, Gamma, delta, params, autocor, aic, bic)
    names(ret) <- c('mllk_optim', 'Gamma', 'delta', 'params', 'autocorrelation', 'AIC', 'BIC')
  } else{
    ret <- list(mod$value, Gamma, delta, params, aic, bic)
    names(ret) <- c('mllk_optim', 'Gamma', 'delta', 'params', 'AIC', 'BIC')
  }
  return(ret)
}

