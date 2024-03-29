% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HMM_sim_functions.R
\name{full_sim_loop}
\alias{full_sim_loop}
\title{Simulation loop}
\usage{
full_sim_loop(
  simulation,
  n_runs,
  dists_fitted,
  p_fitted,
  n_states_fitted,
  n_samples_simulated,
  extract_aic_bic = FALSE,
  ...
)
}
\arguments{
\item{simulation}{Function that generates one simulation.}

\item{n_runs}{Number of models that should be created in the loop.}

\item{dists_fitted}{Vector of distributions of fitted model in R-jargon.}

\item{p_fitted}{Vector of degrees of autoregression of fitted models (one value for every state).}

\item{n_states_fitted}{Number of states of the fitted models.}

\item{n_samples_simulated}{Number of samples simulated in each model.}

\item{extract_aic_bic}{bool, indicates if AIC and BIC of the models should also be saved.}

\item{...}{Input parameters for the simulation function.}
}
\value{
List of model statistics (estimated parameters, accuracies, summary of models).
}
\description{
Simulate a model \code{n_runs} times and refit a certain model to the generated data.
Record the estimated parameters, the accuracies of global decoding and a summary of the models.
}
