# 2022-06-08 
# Functions for density calculation for different distributions

#' Density of gamma distribution
#'
#' Function to calculate the density of different values for a gamma distributed random variable 
#' (with possible autocorrelation in \eqn{mu} and \eqn{sigma}, due to the constant coefficient of variation).
#' 
#' @param x Data.
#' @param theta Named list of parameters of gamma distribution.
#' @param autocor_ind Matrix, indices of to consider when computing the autocorrelated parameters.
#' @param autocor Vector, j-th row of autocor matrix, with one value for each time lag considered.
#' @param p Int, degree of autocorrelation.
#' 
#' @return Densities of gamma distribution.
#' 
#' @export
#' @rdname dens_gamma
dens_gamma <- function(x, theta, autocor_ind=0, autocor=0, p=0){
  cv = theta$sigma/theta$mu
  if (length(theta)==3){ # zero inflation included
    zero_inf = theta$zero_inf
  } else{
    zero_inf = 0
  }
  
  if (p>0){
    # attention: mu and sigma must stay > 0 -> but this is the case by construction
    # this is the autoregression structure
    #mu_auto <- pmax(1e-16, c(#rep(NA,p), # AR(p)
    #             ((1-sum(autocor))*theta$mu + 
    #                as.vector(autocor_ind%*%autocor)))) # matmul of values with autocor coefficient
    mu_auto <- (1-sum(autocor))*theta$mu + 
         as.vector(autocor_ind%*%autocor) # matmul of values with autocor coefficient
    sigma_auto <- cv*mu_auto # calculate sigma using ccv
    
    res = ifelse(x>0, # take zero_inf for values equal to zero
           (1-zero_inf)*dgamma(x, shape = mu_auto^2/sigma_auto^2, scale = sigma_auto^2/mu_auto),
           zero_inf)
    return(res)
  } else{
    res = ifelse(x>0, # take zero_inf for values equal to zero
                 (1-zero_inf)*dgamma(x, shape = theta$mu^2/theta$sigma^2, scale = theta$sigma^2/theta$mu),
                 zero_inf)
    return(res)
  }
}


#' Density of von Mises distribution
#'
#' Function to calculate the density of different values for a von Mises distributed random variable.
#' 
#' @param x Data
#' @param theta Named list of parameters of von Mises distribution.
#' @param autocor_ind Matrix, indices of to consider when computing the autocorrelated parameters.
#' @param autocor Vector, j-th row of autocor matrix, with one value for each time lag considered.
#' @param p Int, degree of autocorrelation.
#' 
#' @return Densities of von Mises distribution.
#' 
#' @export
#' @rdname dens_vm
dens_vm <- function(x, theta, autocor_ind=0, autocor=0, p=0){
  require(CircStats)
  if (p>0){
    ## THIS NEEDS WORK, project [-pi,pi] onto real numbers, to calculate mean properly
    #mu_auto <- c(#rep(NA,p), # AR(p)
    #           ((1-sum(autocor))*theta$mu + 
    #              as.vector(autocor_ind%*%autocor))) # matmul of values with autocor coefficient
    # Alternative for mu_auto:
    mu_auto <- Arg(
      (1-sum(autocor))*exp(1i*theta$mu) + 
        as.vector(exp(1i*autocor_ind)%*%autocor)
      )
    # should work like this
  return(dvm(x, mu_auto, theta$kappa))
  } else{
    return(dvm(x, theta$mu, theta$kappa))
  }
}

#' Density of normal distribution
#'
#' Function to calculate the density of different values for a normal distributed random variable
#' (currently only with possible autocorrelation in parameter \eqn{mu}).
#' 
#' @param x Data
#' @param theta Named list of parameters of normal distribution.
#' @param autocor_ind Matrix, indices of to consider when computing the autocorrelated parameters.
#' @param autocor Vector, j-th row of autocor matrix, with one value for each time lag considered.
#' @param p Int, degree of autocorrelation.
#' 
#' @return Densities of normal distribution.
#' 
#' @export
#' @rdname dens_norm
dens_norm <- function(x, theta, autocor_ind=0, autocor=0, p=0){
  if (p>0){
    mu_auto <- c(#rep(NA,p), # AR(p)
                 ((1-sum(autocor))*theta$mu + 
                    as.vector(autocor_ind%*%autocor))) # matmul of values with autocor coefficient
    return(dnorm(x, mu_auto, theta$sigma))
  } else{
    return(dnorm(x, theta$mu, theta$sigma))
  }
}

