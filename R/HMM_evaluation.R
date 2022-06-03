# 2022-06-03 
# Functions that evaluate different HMMs


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
