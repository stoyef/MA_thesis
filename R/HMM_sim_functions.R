# 2022-05-30
# Functions for simulating and evaluating HMMs with or without autocorrelation



#' Simulate gamma HMMs with and without autocorrelation 
#'
#' Top to bottom wrapper function: Simulate from a specified HMM (autocorrelation 
#' only in parameter \eqn{\mu}), fit a specified
#' HMM to the resulting data and return the output. Optionally: Also return the 
#' decoded states using the Viterbi algorithm. Optionally: Also plot the resulting 
#' density functions of the weighted gamma distributions.
#' 
#' @param model_sim Int, Form of simulated data: Degree of autocorrelation 
#'                       (0 - no autocorrelation).
#' @param model_fit Int, Form of fitted model: Degree of autocorrelation 
#'                       (0 - no autocorrelation).
#' @param N_sim Number of states of the HMM used for data generation.
#' @param N_fit Number of states of the HMM used for data generation.
#' @param n_samples Number of samples generated.
#' @param Gamma_sim Full transition probability matrix of simulated data.
#' @param delta_sim Initial distribution of simulated data.
#' @param mu_sim Parameter vector for mu of simulated data.
#' @param sigma_sim Parameter vector for sigma of simulated data.
#' @param autocor_sim Vector of autocorrelation coefficients of simulated data (if there are any).
#' @param estimate_states Bool, determines if states are estimated and returned
#'                        using Viterbi.
#' @param plot_it Bool, determines if resulting densities are plotted.    
#' 
#' @return List of Fitted model and its parameters (and optional the decoded states).
#' 
#' @export
#' @rdname gamma_simulation_mu     
#' 
gamma_simulation_mu <- function(model_sim, model_fit, N_sim, N_fit, n_samples, 
                             Gamma_sim, delta_sim, mu_sim, sigma_sim, autocor_sim=0,
                             estimate_states=TRUE,plot_it=TRUE){
  
  if (model_sim>0){
    ## prepare autocor for sampling -> matrix
    autocor_matrix <- matrix(autocor_sim, ncol=model_sim, byrow=TRUE) # matrix for easier handling later on
  }
  
  ## simulated data -> contained in simulated_data$data
  if (model_sim==0){
    simulated_data <- sample_hmm_normal(n_samples, delta_sim, Gamma_sim, N_sim,
                                        mu_sim, sigma_sim)
  } else if (model_sim==1){
    simulated_data <- sample_hmm_ar1_mu(n_samples, delta_sim, Gamma_sim, N_sim,
                                           mu_sim, sigma_sim, autocor_matrix)
  } else if (model_sim>1){
    simulated_data <- sample_hmm_arp_mu(n_samples, delta_sim, Gamma_sim, N_sim,
                                     mu_sim, sigma_sim, autocor_matrix,p=model_sim)
  } else{
    return("Wrong input for parameter model_sim.")
  }
  
  # Parameters for model fitting
  theta = c(
    rep(-2,N_fit*(N_fit-1)), # TPM
    rep(0.1,N_fit*model_fit)/(N_fit*model_fit+1), # autocor, n_states*p, regularization to avoid sum > 1
    quantile(simulated_data$data, seq(0,1,length=N_fit)), # mu 
    rep(sd(simulated_data$data),N_fit) # sigma
  )
  
  if (model_fit==0){
    theta.star = c(
      theta[1:(N_fit*(N_fit-1))],
      log(theta[N_fit*(N_fit-1)+1:(2*N_fit)])
    )
  } else{
    theta.star = c(
      theta[1:(N_fit*(N_fit-1))],
      qlogis(theta[(N_fit*(N_fit-1))+1:(N_fit*model_fit)]),
      log(theta[N_fit*(N_fit-1)+(N_fit*model_fit)+1:(2*N_fit)])
    )
  }

  if (model_fit==0){ # no autocorrelation in fitted model
    fitted_model <- fit_arp_model_gamma_mu(mllk_hmm, simulated_data$data,
                                  theta.star, N=N_fit, p=model_fit)
  } else if (model_fit>0){ # there is autocorrelation in fitted model
    fitted_model <- fit_arp_model_gamma_mu(mllk_arp_mu, simulated_data$data,
                                  theta.star, N=N_fit, p=model_fit)
  } else{
    return("Wrong input for parameter model_fit.")
  }
  
  ## Error handling, if optim() in fit_arp_model_gamma_mu didn't work
  # Then we just return an empty output.
  if (anyNA(fitted_model)){
    return(NA)
  }
  
  # Plot, if wanted
  if (plot_it){
    plot_fitted_gamma_dist(simulated_data$data, fitted_model$mu, fitted_model$sigma,
                           fitted_model$delta)
  }
  
  # Viterbi, if wanted
  if (estimate_states){
    estimated_states <- viterbi_arp_mu(simulated_data$data, fitted_model$Gamma,
                                    fitted_model$delta, fitted_model$autocor,
                                    fitted_model$mu, fitted_model$sigma,
                                    N_fit, model_fit)
    
    ret <- list(simulated_data, fitted_model, estimated_states)
    names(ret) <- c('simulated_model','fitted_model', 'viterbi_states')
    return(ret)
  } else{
    ret <- list(simulated_data, fitted_model)
    names(ret) <- c('simulated_model','fitted_model')
    return(ret)
  }
}

