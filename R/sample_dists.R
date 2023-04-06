# 2022-06-10
# Functions for random sample of different distributions

#' Sample from gamma distribution
#'
#' Function to sample from gamma distributed random variable.
#' 
#' @param n Number of samples.
#' @param mu Parameter \eqn{\mu}.
#' @param sigma Parameter \eqn{\sigma}.
#' 
#' @return Sample(s)
#' 
#' @export
#' @rdname sample_gamma
sample_gamma <- function(n, mu, sigma,cv){
  return(rgamma(n, shape=mu^2/sigma^2, scale=sigma^2/mu))
}

#' Sample from von Mises distribution
#'
#' Function to sample from von Mises distributed random variable.
#' 
#' @param n Number of samples.
#' @param mu Parameter \eqn{\mu}.
#' @param kappa Parameter \eqn{\kappa}.
#' 
#' @return Sample(s)
#' 
#' @export
#' @rdname sample_vonMises
#' @import CircStats
sample_vm <- function(n, mu, kappa){
  #require(CircStats)
  return(rvm(n, mean=mu+pi, k=kappa)-pi)
}

#' Sample from normal distribution
#'
#' Function to sample from normal distributed random variable 
#' 
#' @param n Number of samples
#' @param mu Parameter \eqn{\mu}.
#' @param sigma Parameter \eqn{\sigma}.
#' 
#' @return Sample(s)
#' 
#' @export
#' @rdname sample_norm
sample_norm <- function(n, mu, sigma){
  return(rnorm(n, mean=mu, sd=sigma))
}

