# 2022-05-30, updated for general (number of) distributions on 2022-06-14
# Updated for full simulation loop function on 2022-06-24
# Functions for simulating and evaluating HMMs with or without autoregression



#' Simulation of AR(p)-HMMs
#'
#' Top to bottom wrapper function: Simulate from a specified AR(p)-HMM, fit a specified
#' AR(p)-HMM to the resulting data and return the output. Optionally: Also return the 
#' decoded states using the Viterbi algorithm. Optionally: Also plot the resulting 
#' densities of the weighted state-dependent distributions using state-stationary values.
#' 
#' @param model_sim Form of simulated data in a list: First entry is vector containing 
#'                  abbreviated names (in R-jargon) of the distributions 
#'                  the data should be sampled from. Second entry is 
#'                  vector of degree of autoregression for each distribution (one value for each state), 
#'                  0 = no autoregression
#' @param model_fit Form of fitted data in a list: First entry is vector containing 
#'                  abbreviated names (in R-jargon) of the distributions (one value for each state)
#'                  to be considered in the likelihood computation. Second entry is 
#'                  vector of degree of autoregression for each distribution, 
#'                  0 = no autoregression
#' @param N_sim Number of states of the HMM used for data generation.
#' @param N_fit Number of states of the HMM used for fitting the data.
#' @param n_samples Number of samples to be generated.
#' @param Gamma_sim Full transition probability matrix of simulated data.
#' @param delta_sim Initial distribution of simulated data.
#' @param param_sim Parameter vector for parameters of the simulated data. In a vector form,
#'                  and in the order that is customary: e.g. for a gamma HMM c(\eqn{\mu_1,\mu_2,\sigma_1,\sigma_2}).
#'                  The order of distributions has to be the same as in model_sim.
#' @param autocor_sim List of list of parameters for the autocorrelation coefficients. 
#'                    Has to match p, in the order \eqn{\phi_{t-p},\dots,\phi_{t-1}}
#'                    where \eqn{\phi} is the vector of autocorrelation coefficients
#'                   for one specific time lag (one value for each state).
#'                   0, if no autocorrelation. Has to respect the order specified in \code{model_sim}.
#' @param estimate_states Bool, determines if states are estimated and returned
#'                        using Viterbi.
#' @param plot_it Bool, determines if resulting densities are plotted.
#' @param mllk Optional, choice of negative log likelihood function (e.g. penalized).
#' @param lambda Optional, choice of complexity penalty (\eqn{\ge 0}).
#' 
#' @return List of fitted model and its parameters (and optional the decoded states).
#' 
#' @export
#' @rdname ar_simulation     
#' 
ar_simulation <- function(model_sim, model_fit, N_sim, N_fit, n_samples, 
                          Gamma_sim, delta_sim, param_sim, autocor_sim=0,
                          estimate_states=TRUE, plot_it=TRUE, mllk=mllk,
                          lambda=0){
  
  simulated_data <- sample_arp(n_samples=n_samples,
                               delta=delta_sim, 
                               Gamma=Gamma_sim, 
                               N=N_sim, 
                               params=param_sim, 
                               autocor=autocor_sim,
                               p=model_sim[[2]], 
                               dists=model_sim[[1]]
  )
  
  ## starting parameters for model fitting
  params <- c()
  for (dist in 1:length(model_fit[[1]])){
    pars_dist <- starting_params_opt(data=simulated_data$data[,dist],dist=model_fit[[1]][dist],N=N_fit)
    params <- c(params,pars_dist)
  }
  if (any(model_fit[[2]]>0)){
    autocor <- c()
    for (dist in 1:length(model_fit[[1]])){
      ac <- as.numeric(acf(simulated_data$data[,dist], plot=F)$acf[2:(sum(model_fit[[2]][(dist-1)*N_fit+1:N_fit])+1)])
      # Attention: The sum of the autoregression coefficients must not be >1 for a state
      # Otherwise, in the Gamma distribution negative values would be possible for mu
      # which leads to an error
      # Therefore, we standardize ac to avoid this
      ac = ac/(sum(ac)+1)
      autocor <- c(autocor, ac)
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
  fitted_model <- fit_arp_model(mllk=mllk, 
                                data=simulated_data$data, 
                                theta.star=theta.star, 
                                N=N_fit, 
                                p_auto=model_fit[[2]], 
                                dists=model_fit[[1]],
                                lambda=lambda)
  
  
  ## Error handling, if optim() in fit function didn't work
  # Then we just return an empty output.
  if (anyNA(fitted_model)){
    return(NA)
  }
  
  # Plot, if wanted
  if (plot_it){
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






#' Simulation loop
#'
#' Simulate a model \code{n_runs} times and refit a certain model to the generated data.
#' Record the estimated parameters, the accuracies of global decoding and a summary of the models.
#' 
#' @param simulation Function that generates one simulation.
#' @param n_runs Number of models that should be created in the loop.
#' @param dists_fitted Vector of distributions of fitted model in R-jargon.
#' @param p_fitted Vector of degrees of autoregression of fitted models (one value for every state).
#' @param n_states_fitted Number of states of the fitted models.
#' @param n_samples_simulated Number of samples simulated in each model.
#' @param extract_aic_bic bool, indicates if AIC and BIC of the models should also be saved.
#' @param ... Input parameters for the simulation function.
#' 
#' @return List of model statistics (estimated parameters, accuracies, summary of models).
#' 
#' @export
#' @rdname full_sim_loop     
#' 
full_sim_loop <- function(simulation, n_runs, dists_fitted, p_fitted, 
                          n_states_fitted, n_samples_simulated, extract_aic_bic=FALSE,...){
  
  # Inputs for simulation function
  args=list(...)
  
  # Create matrices where the estimated parameters are stored
  # Currently only works for distributions with 2 parameters
  for (dist in 1:length(dists_fitted)){
    assign(paste("estimated_",dist,"_param_1", sep=""), 
           matrix(NA,nrow=n_runs,ncol=n_states_fitted))
    assign(paste("estimated_",dist,"_param_2", sep=""), 
           matrix(NA,nrow=n_runs,ncol=n_states_fitted))
    for (state in 1:n_states_fitted){
      assign(paste("estimated_",dist,"state_",state,"_autocor", sep=""),
             matrix(NA,nrow=n_runs,ncol=p_fitted[(dist-1)*n_states_fitted+state]))
    }
  }
  true_states = matrix(NA,nrow=n_runs,ncol=n_samples_simulated)
  estimated_states = matrix(NA,nrow=n_runs,ncol=n_samples_simulated)
  aics = rep(NA,n_runs)
  bics = rep(NA,n_runs)
  logLike = rep(NA,n_runs)
  
  start_time = Sys.time()
  
  ## Update 2023-02-02 -- use of optimParallel instead of optim makes parallelization
  ## of simulation runs unnecessary (all cores are working anyway)
  n_its = 0 # counter of successful iterations
  results = list()
  
  while(n_its < n_runs){ # run as long as all models are fitted
    current_time=Sys.time()
    sim = do.call(ar_simulation, args)
    if (length(sim)>0){ # if sim was sucessful
      n_its = n_its + 1
      # insert estimated parameters in matrices by only accessing string values
      # -> weird workaround with get() and temporary matrix
      for (dist in 1:length(dists_fitted)){
        # 1st parameter
        h = get(paste("estimated_",dist,"_param_1", sep=""))
        h[n_its,] = sim$fitted_model$params[[dist]][[1]]
        assign(paste("estimated_",dist,"_param_1", sep=""), h)
        # 2nd parameter
        h = get(paste("estimated_",dist,"_param_2", sep=""))
        h[n_its,] = sim$fitted_model$params[[dist]][[2]]
        assign(paste("estimated_",dist,"_param_2", sep=""), h)
        # autocorrelation
        for (state in 1:n_states_fitted){
          h = get(paste("estimated_",dist,"state_",state,"_autocor", sep=""))
          h[n_its,] = sim$fitted_model$autocorrelation[[dist]][[state]]
          assign(paste("estimated_",dist,"state_",state,"_autocor", sep=""), h)
        }
      }
      true_states[n_its,] = sim$simulated_model$states
      estimated_states[n_its,] = sim$viterbi_states
      if (extract_aic_bic){
        # AIC
        aics[n_its] = sim$fitted_model$AIC
        # BIC
        bics[n_its] = sim$fitted_model$BIC
        # log-likelihood
        logLike[n_its] = - sim$fitted_model$mllk_optim
      }
    }
  }
  
  
  ### OLD:
  ###
  ## separate code for single and multicore computations
  #if (multicore){
  #  require(parallel)
  #  n_cores=detectCores()
  #  # we need a wrapper function for mclapply
  #  sim_wrap <- function(iteration){
  #    current_time=Sys.time()
  #    sim <- do.call(ar_simulation,args)
  #    return(sim)
  #  }
  #  
  #  n_its = 0 # counter of successful iterations
  #  while(length(which(is.na(estimated_1_param_1[,1])))>0){ # run as long as all models are fitted
  #    iterations = 1:length(which(is.na(estimated_1_param_1[,1]))) # missing iterations
  #    # parallelization with mclapply:
  #    results <- mclapply(iterations,
  #                        sim_wrap,
  #                        mc.cores = n_cores
  #    )
  #    
  #    # check which iterations were successful
  #    successful_iterations_this_time = unlist(lapply(results, function(it) length(it)>1))
  #    successful_results_this_time = results[successful_iterations_this_time]
  #    
  #    if (sum(successful_iterations_this_time)>0){ # execute loop only, if some fits were successful
  #      for (i in n_its+1:sum(successful_iterations_this_time)){ # fill up the matrices top to bottom
  #        # insert estimated parameters in matrices by only accessing string values
  #        # -> weird workaround with get() and temporary matrix
  #        for (dist in 1:length(dists_fitted)){
  #          # 1st parameter
  #          h = get(paste("estimated_",dist,"_param_1", sep=""))
  #          h[i,] = successful_results_this_time[[i-n_its]]$fitted_model$params[[dist]][[1]]
  #          assign(paste("estimated_",dist,"_param_1", sep=""), h)
  #          # 2nd parameter
  #          h = get(paste("estimated_",dist,"_param_2", sep=""))
  #          h[i,] = successful_results_this_time[[i-n_its]]$fitted_model$params[[dist]][[2]]
  #          assign(paste("estimated_",dist,"_param_2", sep=""), h)
  #          # autocorrelation
  #          for (state in 1:n_states_fitted){
  #            h = get(paste("estimated_",dist,"state_",state,"_autocor", sep=""))
  #            h[i,] = successful_results_this_time[[i-n_its]]$fitted_model$autocorrelation[[dist]][[state]]
  #            assign(paste("estimated_",dist,"state_",state,"_autocor", sep=""), h)
  #          }
  #        }
  #        true_states[i,] = successful_results_this_time[[i-n_its]]$simulated_model$states
  #        estimated_states[i,] = successful_results_this_time[[i-n_its]]$viterbi_states
  #        if (extract_aic_bic){
  #          # AIC
  #          aics[i] = successful_results_this_time[[i-n_its]]$fitted_model$AIC
  #          # BIC
  #          bics[i] = successful_results_this_time[[i-n_its]]$fitted_model$BIC
  #          # log-likelihood
  #          logLike[i] = - successful_results_this_time[[i-n_its]]$fitted_model$mllk_optim
  #        }
  #      }
  #      n_its = n_its+sum(successful_iterations_this_time)
  #    }
  #  }
  #  
  #} else{ # single core
  #  while(length(which(is.na(estimated_1_param_1[,1])))>0){ # run as long as all models are fitted
  #    for (i in which(is.na(estimated_1_param_1[,1]))){ # re-run only models that failed last time
  #      current_time = Sys.time()
  #      sim <- do.call(ar_simulation,args)
  #      cat(i,'/',n_runs,' (', Sys.time()-current_time,') \n',sep="")
  #      
  #      # error handling, skip iteration if optim() in fit function didn't work
  #      if(anyNA(sim)){
  #        next
  #      }
  #      
  #      # insert estimated parameters in matrices by only accessing string values
  #      # -> weird workaround with get() and temporary matrix
  #      for (dist in 1:length(dists_fitted)){
  #        # 1st parameter
  #        h = get(paste("estimated_",dist,"_param_1", sep=""))
  #        h[i,] = sim$fitted_model$params[[dist]][[1]]
  #        assign(paste("estimated_",dist,"_param_1", sep=""), h)
  #        # 2nd parameter
  #        h = get(paste("estimated_",dist,"_param_2", sep=""))
  #        h[i,] = sim$fitted_model$params[[dist]][[2]]
  #        assign(paste("estimated_",dist,"_param_2", sep=""), h)
  #        # autocorrelation
  #        for (state in 1:n_states_fitted){
  #          h = get(paste("estimated_",dist,"state_",state,"_autocor", sep=""))
  #          h[i,] = sim$fitted_model$autocorrelation[[dist]][[state]]
  #          assign(paste("estimated_",dist,"state_",state,"_autocor", sep=""), h)
  #        }
  #      }
  #      true_states[i,] = sim$simulated_model$states
  #      estimated_states[i,] = sim$viterbi_states
  #      if (extract_aic_bic){
  #        # AIC
  #        aics[i] = sim$fitted_model$AIC
  #        # BIC
  #        bics[i] = sim$fitted_model$BIC
  #        # log-likelihood
  #        logLike[i] = - sim$fitted_model$mllk_optim
  #      }
  #    }
  #  }
  #}
  
  elapsed_time = Sys.time()-start_time
  cat("Total simulation time:", elapsed_time, "\n")
  
  # accuracies
  acc = rep(NA,n_runs)
  for (i in (1:n_runs)){
    acc[i]=sum(true_states[i,] == estimated_states[i,])/n_samples_simulated
  }
  
  param_estimates <- list()
  autocor_estimates <- list()
  for (dist in 1:length(dists_fitted)){
    param_estimates[[paste("estimated_",dist,"_param_1", sep="")]] = get(paste("estimated_",dist,"_param_1", sep=""))
    param_estimates[[paste("estimated_",dist,"_param_2", sep="")]] = get(paste("estimated_",dist,"_param_2", sep=""))
    for (state in 1:n_states_fitted){
      autocor_estimates[[paste("estimated_",dist,"state_",state,"_autocor", sep="")]] = get(paste("estimated_",dist,"state_",state,"_autocor", sep=""))
    }
  }
  
  if (extract_aic_bic){
    ret = list(param_estimates, autocor_estimates, acc, logLike, aics, bics)
    names(ret) = c("estimated_parameters", "estimated_autocorrelation", "decoding_accuracies", "log-Likelihood", "AIC", "BIC")
  } else{
    ret = list(param_estimates, autocor_estimates, acc)
    names(ret) = c("estimated_parameters", "estimated_autocorrelation", "decoding_accuracies")
  }
  return(ret)
}

