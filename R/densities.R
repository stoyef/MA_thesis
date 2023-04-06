# 2022-06-08 
# Functions for density calculation for different distributions

#' Density of gamma distribution
#'
#' Function to calculate the density of different values for a gamma distributed random variable 
#' (with possible autoregression in \eqn{mu} and \eqn{sigma}, due to the constant coefficient of variation).
#' 
#' @param x Data.
#' @param theta Named list of parameters of gamma distribution.
#' @param autocor_ind Matrix, indices to consider when computing the autoregressive parameters.
#' @param autocor Vector, j-th row of autocor matrix, with one value for each time lag considered.
#' @param p Int, degree of autoregression.
#' 
#' @return Vector of densities of gamma distribution.
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
#' Function to calculate the density of different values for a von Mises distributed random variable
#' (with possible autoregression in \eqn{mu}).
#' 
#' @param x Data.
#' @param theta Named list of parameters of von Mises distribution.
#' @param autocor_ind Matrix, indices to consider when computing the autoregressive parameters.
#' @param autocor Vector, j-th row of autocor matrix, with one value for each time lag considered.
#' @param p Int, degree of autoregression
#' 
#' @return Vector of densities of von Mises distribution.
#' 
#' @export
#' @rdname dens_vm
#' @import CircStats
dens_vm <- function(x, theta, autocor_ind=0, autocor=0, p=0){
  #require(CircStats)
  if (p>0){
    mu_auto <- Arg(
      (1-sum(autocor))*exp(1i*theta$mu) + 
        as.vector(exp(1i*autocor_ind)%*%autocor)
      )
  return(dvm(x, mu_auto, theta$kappa))
  } else{
    return(dvm(x, theta$mu, theta$kappa))
  }
}

#' Density of normal distribution
#'
#' Function to calculate the density of different values for a normal distributed random variable
#' (with possible autoregression in parameter \eqn{mu}).
#' 
#' @param x Data.
#' @param theta Named list of parameters of normal distribution.
#' @param autocor_ind Matrix, indices to consider when computing the autoregressive parameters.
#' @param autocor Vector, j-th row of autocor matrix, with one value for each time lag considered.
#' @param p Int, degree of autoregression
#' 
#' @return Vector of densities of normal distribution.
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

