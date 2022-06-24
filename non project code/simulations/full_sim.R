## 2022-06-24 Full simulation loop function

library(MasterThesis)

par(mfrow=c(1,2))


# sim2, fit2
full_sim_1 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 5,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(2,2),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(2,2)),
  model_fit=list(c('gamma','vm'),c(2,2)),
  N_sim=2,
  N_fit=2,
  n_samples=2000,
  Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),ncol=2),
  delta_sim = c(0.5,0.5),
  param_sim = c(20,40,5,7,
                0,0,2,12),
  autocor_sim = list(matrix(c(0.15,0.3,0.15,0.4),ncol=2,byrow=TRUE),
                     matrix(c(0.2,0.3,0.2,0.4),ncol=2,byrow=TRUE)),
  estimate_states = TRUE,
  plot_it = TRUE
)

full_sim_1


# sim0, fit2
full_sim_2 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 5,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(2,2),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(0,0)),
  model_fit=list(c('gamma','vm'),c(2,2)),
  N_sim=2,
  N_fit=2,
  n_samples=2000,
  Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),ncol=2),
  delta_sim = c(0.5,0.5),
  param_sim = c(20,40,5,7,
                0,0,2,12),
  autocor_sim = 0,#list(matrix(c(0.15,0.3,0.15,0.4),ncol=2,byrow=TRUE),
                 #     matrix(c(0.2,0.3,0.2,0.4),ncol=2,byrow=TRUE)),
  estimate_states = TRUE,
  plot_it = TRUE
)

full_sim_2


# sim2, fit0
full_sim_3 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 5,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(0,0),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(2,2)),
  model_fit=list(c('gamma','vm'),c(0,0)),
  N_sim=2,
  N_fit=2,
  n_samples=2000,
  Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),ncol=2),
  delta_sim = c(0.5,0.5),
  param_sim = c(20,40,5,7,
                0,0,2,12),
  autocor_sim = list(matrix(c(0.15,0.3,0.15,0.4),ncol=2,byrow=TRUE),
                     matrix(c(0.2,0.3,0.2,0.4),ncol=2,byrow=TRUE)),
  estimate_states = TRUE,
  plot_it = TRUE
)

full_sim_3
