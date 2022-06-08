# 2022-06-08
# Function to construct nice HMM parameters in named-list form

#' Construct parameter object for HMM likelihood evaluation
#'
#' Construct a nice named list that contains all relevant parameters of an HMM.
#' This gets used in the Log-Likelihood computation to generalize to arbitrary 
#' distributions.
#' 
#' @param Gamma Full TPM.
#' @param autocor Vector or matrix of autocorrelation coefficients in suitable format.
#' @param ... List of other parameters, specified in "param" (e.g. mu, sigma, depending on distribution). 
#'            Here, notational conventions have to be considered when supplying the distributional parameters
#' 
#' @return Named list of all parameters
#' 
#' @export
#' @rdname nice_params
nice_params <- function(Gamma, autocor, ...){
  args = list(Gamma=Gamma,autocor=autocor,...)
  return(args)
}