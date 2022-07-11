### 2022-07-11
### Multicore computation using package parallel
### and function mclapply

library(parallel)

n_cores = detectCores()
n_cores
## My machine:
## Apple M1 Max SOC
## 10 Core CPU

## For each model we need a wrapper functions, otherwise the function
## mclapply does not work
ar_simulation_wrap_11 <- function(iteration){
  res <- ar_simulation(
    model_sim=list(c('gamma','vm'),c(1,1)),
    model_fit=list(c('gamma','vm'),c(1,1)),
    N_sim=2,
    N_fit=2,
    n_samples=2000,
    Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),ncol=2),
    delta_sim = c(0.5,0.5),
    param_sim = c(20,40,5,7,
                  0,0,2,12),
    autocor_sim = list(matrix(c(0.3,0.7),ncol=1),
                       matrix(c(0.4,0.7),ncol=1)),
    estimate_states = TRUE,
    plot_it = FALSE # plots don't work with parallelization
  )
  return(res)
}

iterations = 3:6
results <- mclapply(iterations,
                    ar_simulation_wrap_11,
                    mc.cores = n_cores
                    )
# works


## let's try the whole deal

# sim0, fit0
full_sim_00 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(0,0),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(0,0)),
  model_fit=list(c('gamma','vm'),c(0,0)),
  N_sim=2,
  N_fit=2,
  n_samples=2000,
  Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),ncol=2),
  delta_sim = c(0.5,0.5),
  param_sim = c(20,40,5,7,
                0,0,2,12),
  autocor_sim = 0,#list(matrix(c(0.15,0.3,0.15,0.4),ncol=2,byrow=TRUE),
  #   matrix(c(0.2,0.3,0.2,0.4),ncol=2,byrow=TRUE)),
  estimate_states = TRUE,
  plot_it = FALSE,
  multicore = TRUE
)







