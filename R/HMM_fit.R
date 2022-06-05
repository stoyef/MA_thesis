# 2022-06-03 
# Functions that fit different HMMs to data


#' Fit an AR(p) HMM to data
#' 
#' Fit an AR(p) HMM with autocorrelation to data 
#' using a specified function to compute the negative 
#' Log-Likelihood. This function gets minimized using the function \code{optim}.
#' It returns the estimated parameters of the fitted model. The distribution of the HMM
#' can be specified using the parameter \code{dist}.
#' 
#' @param mllk Negative log-Likelihood function that should be minimized.
#' @param data Data that should be fitted using the HMM.
#' @param theta.star Unconstrained parameter vector (has to be provided in suitable form).
#' @param N Number of states.
#' @param p Degree of autocorrelation, 0 is possible for no autocorrelation.
#' @param dist Kind of distribution in the HMM. One of ['gamma', 'von Mises'].
#' 
#' @return List, containing minimal value of negative log-Likelihood, Gamma, 
#'         delta, (autocorrelation, depending on degree), mu, sigma.
#' 
#' @export
#' @rdname fit_arp_model
fit_arp_model <- function(mllk, data, theta.star, N, p, dist){
  # Error handling: sometimes the optim function does not work, because of some
  # error. We want to be notified that there is an error, but the execution should 
  # not be interrupted (important for for loops that use this function)
  skip = FALSE
    tryCatch(
      mod <- optim(par=theta.star, fn=mllk, method='L-BFGS-B',
                   N=N,p=p,x=data),
      error=function(e){cat("ERROR: optim() failed. Continue to next iteration.\n")
        skip <<- TRUE
      }
      
    )
  
  # leave function, if optim() didn't work
  if (skip){
    return(NA)
  }
  
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(mod$par[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  if (p>0){
    if (dist=='gamma'){
      # mu and sigma 
      mu <- exp(mod$par[(N-1)*N+1:N])
      sigma <- exp(mod$par[(N-1)*N+N+1:N])
      # autocorrelation
      autocor <- plogis(mod$par[(N-1)*N+2*N+1:(p*N)])
      ret <- list(mod$value, Gamma, delta, autocor, mu, sigma)
      names(ret) <- c('mllk', 'Gamma', 'delta', 'autocorrelation', 'mu', 'sigma')
      return(ret)
    }
    if (dist=='von Mises'){
      # mu and kappa 
      mu <- Arg(mod$par[(N-1)*N+1:N]+1i*mod$par[(N-1)*N+N+1:N])
      kappa <- sqrt(mod$par[(N-1)*N+1:N]^2+mod$par[(N-1)*N+N+1:N]^2)
      # autocorrelation
      autocor <- plogis(mod$par[(N-1)*N+2*N+1:(p*N)])
      ret <- list(mod$value, Gamma, delta, autocor, mu, kappa)
      names(ret) <- c('mllk', 'Gamma', 'delta', 'autocorrelation', 'mu', 'kappa')
      return(ret)
    }
  } else{
    if (dist=='gamma'){
      # mu and sigma 
      mu <- exp(mod$par[(N-1)*N+1:N])
      sigma <- exp(mod$par[(N-1)*N+N+1:N])
      ret <- list(mod$value, Gamma, delta, mu, sigma)
      names(ret) <- c('mllk', 'Gamma', 'delta', 'mu', 'sigma')
      return(ret)
    }
    if (dist=='von Mises'){
      # mu and kappa 
      mu <- Arg(mod$par[(N-1)*N+1:N]+1i*mod$par[(N-1)*N+N+1:N])
      kappa <- sqrt(mod$par[(N-1)*N+1:N]^2+mod$par[(N-1)*N+N+1:N]^2)
      ret <- list(mod$value, Gamma, delta, mu, kappa)
      names(ret) <- c('mllk', 'Gamma', 'delta', 'mu', 'kappa')
      return(ret)
    }
  }
}
