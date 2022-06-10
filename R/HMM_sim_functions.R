# 2022-05-30
# Functions for simulating and evaluating HMMs with or without autocorrelation



#' Simulate HMMs with and without autocorrelation 
#'
#' Top to bottom wrapper function: Simulate from a specified HMM, fit a specified
#' HMM to the resulting data and return the output. Optionally: Also return the 
#' decoded states using the Viterbi algorithm. Optionally: Also plot the resulting 
#' density functions of the weighted distributions.
#' 
#' @param model_sim Form of simulated data in a vector: First entry is vector containing 
#'                  abbreviated names (in R-jargon) of the distributions 
#'                  the data should be sampled from. Second entry is 
#'                  vector of degree of autocorrelation for each distribution, 
#'                  0=no autocorrelation.
#' @param model_fit Form of fitted data in a vector: First entry is vector containing 
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
#' @param autocor_sim Vector of autocorrelation coefficients of simulated data (if there are any,
#'                    one value for each time lag and state). The order of distributions has
#'                    to be the same as in model_sim.
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
  
  if (as.integer(model_sim[2])>0){
    ## prepare autocor for sampling -> matrix
    autocor_matrix <- matrix(autocor_sim, ncol=as.integer(model_sim[2]), byrow=TRUE) # matrix for easier handling later on
  }
  
  ## simulated data -> contained in simulated_data$data
  if (model_sim[1] == 'gamma'){
    simulated_data <- sample_gamma_arp(n_samples, delta_sim, Gamma_sim, N_sim,
                                       param_sim[1:N_sim], param_sim[(N_sim+1):(2*N_sim)], 
                                       autocor_matrix, p=as.integer(model_sim[2]))
  } else if (model_sim[1] == 'von Mises'){
    simulated_data <- sample_vonMises_arp(n_samples, delta_sim, Gamma_sim, N_sim,
                                          param_sim[1:N_sim], param_sim[(N_sim+1):(2*N_sim)],
                                          autocor_matrix, p=as.integer(model_sim[2]))
  } else{
    return("Wrong input for parameter model_sim.")
  }
  
  # Parameters for model fitting, differenciate for the distributions
  if (model_fit[1]=='gamma'){
    theta = c(
      rep(-2,N_fit*(N_fit-1)), # TPM
      quantile(simulated_data$data, seq(0.1,0.9,length=N_fit)), # mu 
      rep(sd(simulated_data$data),N_fit), # sigma
      rep(0.1,N_fit*as.integer(model_fit[2]))/(N_fit*as.integer(model_fit[2])+1) # autocor, n_states*p, regularization to avoid sum > 1
    )
    
    if (as.integer(model_fit[2])==0){
      theta.star = c(
        theta[1:(N_fit*(N_fit-1))],
        log(theta[N_fit*(N_fit-1)+1:(2*N_fit)])
      )
    } else{
      theta.star = c(
        theta[1:(N_fit*(N_fit-1))],
        log(theta[N_fit*(N_fit-1)+1:(2*N_fit)]),
        qlogis(theta[(N_fit*(N_fit-1))+(2*N_fit)+1:(N_fit*as.integer(model_fit[2]))])
      )
    }
    
    # fit the model
    fitted_model <- fit_arp_model(mllk_gamma_arp, simulated_data$data,
                                  theta.star, N=N_fit, p=as.integer(model_fit[2]), dist=model_fit[1])
    
  } else if (model_fit[1]=='von Mises'){
    theta = c(
      rep(-2,N_fit*(N_fit-1)), # TPM
      quantile(simulated_data$data, seq(0.25,0.75,length=N_fit)), # mu 
      rep(est.kappa(simulated_data$data),N_fit), # kappa
      rep(0.1,N_fit*as.integer(model_fit[2]))/(N_fit*as.integer(model_fit[2])+1) # autocor, n_states*p, regularization to avoid sum > 1
    )
    
    if (model_fit[2]==0){
      theta.star = c(
        theta[1:(N_fit*(N_fit-1))],
        theta[N_fit*(N_fit-1)+1:N_fit] * cos(theta[N_fit*(N_fit-1)+N_fit+1:N_fit]), # mean
        theta[N_fit*(N_fit-1)+1:N_fit] * sin(theta[N_fit*(N_fit-1)+N_fit+1:N_fit]), # kappa
      )
    } else{
      theta.star = c(
        theta[1:(N_fit*(N_fit-1))],
        theta[N_fit*(N_fit-1)+1:N_fit] * cos(theta[N_fit*(N_fit-1)+N_fit+1:N_fit]), # mean
        theta[N_fit*(N_fit-1)+1:N_fit] * sin(theta[N_fit*(N_fit-1)+N_fit+1:N_fit]), # kappa
        qlogis(theta[(N_fit*(N_fit-1))+(2*N_fit)+1:(N_fit*as.integer(model_fit[2]))])
      )
    }
    
    # fit the model
    fitted_model <- fit_arp_model(mllk_vonMises_arp, simulated_data$data,
                                  theta.star, N=N_fit, p=as.integer(as.integer(model_fit[2])), dist=model_fit[1])

  }
  
  ## Error handling, if optim() in fit function didn't work
  # Then we just return an empty output.
  if (anyNA(fitted_model)){
    return(NA)
  }
  
  # Plot, if wanted
  if (plot_it){
    if (model_fit[1]=='gamma'){
      param <- c(fitted_model$mu,fitted_model$sigma)
    } else if (model_fit[1]=='von Mises'){
      param <- c(fitted_model$mu,fitted_model$kappa)
    }
    plot_fitted_dist(simulated_data$data, model_fit[1], param, N_fit,
                           fitted_model$delta)
    }
  
  # Viterbi, if wanted
  if (estimate_states){
    if (model_fit[1]=='gamma'){
      estimated_states <- viterbi_gamma_arp(simulated_data$data, fitted_model$Gamma,
                                      fitted_model$delta, fitted_model$autocor,
                                      fitted_model$mu, fitted_model$sigma,
                                      N_fit, as.integer(model_fit[2]))
      
      ret <- list(simulated_data, fitted_model, estimated_states)
      names(ret) <- c('simulated_model','fitted_model', 'viterbi_states')
      return(ret)
    } else if (model_fit[1]=='von Mises'){
      estimated_states <- viterbi_vonMises_arp(simulated_data$data, fitted_model$Gamma,
                                            fitted_model$delta, fitted_model$autocor,
                                            fitted_model$mu, fitted_model$kappa,
                                            N_fit, as.integer(model_fit[2]))
      
      ret <- list(simulated_data, fitted_model, estimated_states)
      names(ret) <- c('simulated_model','fitted_model', 'viterbi_states')
      return(ret)
      }
    } else{
      ret <- list(simulated_data, fitted_model)
      names(ret) <- c('simulated_model','fitted_model')
      return(ret)
    }
}

