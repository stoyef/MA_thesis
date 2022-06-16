# 2022-05-30, updated for general (number of) distributions on 2022-06-14
# Functions for simulating and evaluating HMMs with or without autocorrelation



#' Simulate HMMs with and without autocorrelation 
#'
#' Top to bottom wrapper function: Simulate from a specified HMM, fit a specified
#' HMM to the resulting data and return the output. Optionally: Also return the 
#' decoded states using the Viterbi algorithm. Optionally: Also plot the resulting 
#' density functions of the weighted distributions.
#' 
#' @param model_sim Form of simulated data in a list: First entry is vector containing 
#'                  abbreviated names (in R-jargon) of the distributions 
#'                  the data should be sampled from. Second entry is 
#'                  vector of degree of autocorrelation for each distribution, 
#'                  0=no autocorrelation.
#' @param model_fit Form of fitted data in a list: First entry is vector containing 
#'                  abbreviated names (in R-jargon) of the distributions 
#'                  to be considered in the Likelihood computation. Second entry is 
#'                  vector of degree of autocorrelation for each distribution, 
#'                  0=no autocorrelation.
#' @param N_sim Number of states of the HMM used for data generation.
#' @param N_fit Number of states of the HMM used for fitting the data.
#' @param n_samples Number of samples to be generated.
#' @param Gamma_sim Full transition probability matrix of simulated data.
#' @param delta_sim Initial distribution of simulated data.
#' @param param_sim Parameter vector for parameters of the simulated data. In a vector form,
#'                  and in the order that is customary: e.g. for a gamma HMM c(\eqn{\mu1,\mu2,\sigma1,\sigma2}).
#'                  The order of distributions has to be the same as in model_sim.
#' @param autocor_sim List of parameter matrix (Dimensions for each matrix: \eqn{n\times p}) for the autocorrelation coefficients. 
#'                    Has to match p, in the order \eqn{\phi_{t-p},\dots,\phi_{t-1}}
#'                    where \eqn{\phi} is the vector of autocorrelation coefficients
#'                   for one specific time lag (one value for each state).
#'                   0, if no autocorrelation. Has to respect the order specified in \code{model_sim}.
#' @param estimate_states Bool, determines if states are estimated and returned
#'                        using Viterbi.
#' @param plot_it Bool, determines if resulting densities are plotted.    
#' 
#' @return List of Fitted model and its parameters (and optional the decoded states).
#' 
#' @export
#' @rdname ar_simulation     
#' 
ar_simulation <- function(model_sim, model_fit, N_sim, N_fit, n_samples, 
                             Gamma_sim, delta_sim, param_sim, autocor_sim=0,
                             estimate_states=TRUE,plot_it=TRUE){
  
  #if (as.integer(model_sim[[2]])>0){
  #  ## prepare autocor for sampling -> matrix
  #  autocor_matrix <- matrix(autocor_sim, ncol=as.integer(model_sim[2]), byrow=TRUE) # matrix for easier handling later on
  #}
  
  simulated_data <- sample_arp(n_samples=n_samples,
                               delta=delta_sim, 
                               Gamma=Gamma_sim, 
                               N=N_sim, 
                               params=param_sim, 
                               autocor=autocor_sim,
                               p=model_sim[[2]], 
                               dists=model_sim[[1]]
                                 )
  ## simulated data -> contained in simulated_data$data

  ## starting parameters for model fitting
  params <- c()
  for (dist in 1:length(model_fit[[1]])){
    pars_dist <- starting_params_opt(data=simulated_data$data[,dist],dist=model_fit[[1]][dist],N=N_fit)
    params <- c(params,pars_dist)
  }
  if (any(model_fit[[2]]>0)){
    autocor <- c()
    for (dist in 1:length(model_fit[[1]])){
      ac <- as.numeric(acf(simulated_data$data[,dist], plot=F)$acf[2:(model_fit[[2]][dist]+1)])
      autocor <- c(autocor, rep(ac,N_fit))
    }
    theta <- c(
      rep(-2,N_fit*(N_fit-1)), #TPM
      params, # dist parameters
      autocor # autocor parameters
    )
  } else{
    theta <- c(
      rep(-2,N_fit*(N_fit-1)), #TPM
      params # dist parameters
    )
  }
  
  theta.star <- starize(theta=theta, N=N_fit, p=model_fit[[2]], dists=model_fit[[1]])
  
  # fit the model
  fitted_model <- fit_arp_model(mllk, 
                    simulated_data$data, 
                    theta.star, 
                    N=N_fit, 
                    p=model_fit[[2]], 
                    dists=model_fit[[1]])
  
  
  ## Error handling, if optim() in fit function didn't work
  # Then we just return an empty output.
  if (anyNA(fitted_model)){
    return(NA)
  }
  
  # Plot, if wanted
  if (plot_it){
    #cat("Generic plot function not implemented yet.")
    for (dist in 1:length(model_fit[[1]])){
      if (model_fit[[1]][dist] == 'norm'){
        cat("Plot for normal distribution not implemented yet.\n")
      }
      param <- fitted_model$params[[dist]]
      plot_fitted_dist(simulated_data$data[,dist], model_fit[[1]][dist], param, N_fit,
                        fitted_model$delta)
      
    }
  }
  
  
  
  
  # Viterbi, if wanted
  if (estimate_states){
    #cat("Generic decoding function not implemented yet.")
    estimated_states <- viterbi_arp(x=simulated_data$data, 
                                    Gamma=fitted_model$Gamma,
                                    delta=fitted_model$delta, 
                                    dists=model_fit[[1]],
                                    autocor=fitted_model$autocor,
                                    params=fitted_model$params,
                                    N=N_fit, 
                                    p=model_fit[[2]])
    
    ret <- list(simulated_data, fitted_model, estimated_states)
    names(ret) <- c('simulated_model','fitted_model', 'viterbi_states')
    return(ret)
    } else{
      ret <- list(simulated_data, fitted_model)
      names(ret) <- c('simulated_model','fitted_model')
      return(ret)
    }
}

