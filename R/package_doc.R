#' MasterThesis: To give all this code some structure
#'
#' This package includes functions that employ the user to create simulation
#' studies with HMMs that contain an AR(p)-structure in the state dependent process.
#' The functions are primarily used in my Master Thesis (hence the package name :D).
#' Of course, everyone is free to use them in any other context.
#' 
#' @section Important functions:
#' \code{\link{sample_arp}} Simulate data from an HMM with AR(p) structure
#' 
#' \code{\link{mllk}} Compute negative log Likelihood of an AR(p) HMM
#' 
#' \code{\link{fit_arp_model}} Fit an AR(p) HMM to data
#' 
#' \code{\link{ar_simulation}} Simulate data from an AR(p) HMM and refit specified model to it
#' 
#' \code{\link{full_sim_loop}} Simulate and refit models iteratively and record model specifications
#'  
#' @references 
#' Zucchini, Walter, Iain L. MacDonald and Roland Langrock (2016). Hidden Markov
#' Models for Time Series: An Introduction Using R. 2nd ed. Boca Raton: Chapman
#' and Hall/CRC. doi: 10.1201/9781420010893.
#'
#'
#' @docType package
#' @name MasterThesis
NULL
#> NULL