#' MasterThesis: To give all this code some structure
#'
#' This package includes functions that employ the user to create simulation
#' studies with AR(p)-HMMs (HMMs with within-state autoregression) and to fit AR(p)-HMMs to real data. 
#' Until now, the motivation of AR(p)-HMMs has been animal movement data, using gamma distributed
#' step lengths and von Mises distributed turning angles. 
#' Therefore, funcitonalities can only be guaranteed to work with bivariate AR(p)-HMMs with one gamma distributed
#' and one von Mises distributed variable (although in principle this can be changed).
#' The functions are primarily used in my Master Thesis (hence the package name :D).
#' Of course, everyone is free to use them in any other context.
#' 
#' @section Important functions:
#' \code{\link{sample_arp}} Simulate data from an AR(p)-HMM 
#' 
#' \code{\link{mllk}} Compute negative log-likelihood of an AR(p)-HMM
#' 
#' \code{\link{fit_arp_model}} Fit an AR(p)-HMM to data
#' 
#' \code{\link{ar_simulation}} Simulate data from an AR(p)-HMM and refit specified AR(p)-HMM to it
#' 
#' \code{\link{full_sim_loop}} Simulate and refit models iteratively and record model specifications
#'  
#' @references 
#' Zucchini, Walter, Iain L. MacDonald and Roland Langrock (2016). Hidden Markov
#' Models for Time Series: An Introduction Using R. 2nd ed. Boca Raton: Chapman
#' and Hall/CRC. doi: 10.1201/9781420010893.
#' 
#' Michelot, Theo, Roland Langrock and Toby A. Patterson (2016). “moveHMM: 
#' an R package for the statistical modelling of animal movement data using 
#' hidden Markov models”. In: Methods in Ecology and Evolution 7.11, pp. 1308–1315.
#'
#'
#' @docType package
#' @name MasterThesis
NULL
#> NULL