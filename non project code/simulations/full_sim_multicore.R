### 2022-07-11
### Multicore computation using package parallel
### and function mclapply

library(parallel)

n_cores = detectCores()
n_cores
## My machine:
## Apple M1 Max SOC, 32 GB shared LPDDR5 memory
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

################################################################################

par(mfrow=c(1,2))

total_begin = Sys.time()
# delete later: stores computation times
comps = rep(NA,16)
c=1



# sim0, fit0
t=Sys.time()
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

comps[c]=Sys.time()-t
c=c+1

# full_sim_00

################################################################################

# sim0, fit1
t=Sys.time()
full_sim_01 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(1,1),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(0,0)),
  model_fit=list(c('gamma','vm'),c(1,1)),
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

comps[c]=Sys.time()-t
c=c+1

#full_sim_01

################################################################################

# sim0, fit2
t=Sys.time()
full_sim_02 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
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
  #   matrix(c(0.2,0.3,0.2,0.4),ncol=2,byrow=TRUE)),
  estimate_states = TRUE,
  plot_it = FALSE,
  multicore = TRUE
)

comps[c]=Sys.time()-t
c=c+1

#full_sim_02

################################################################################

# sim0, fit3
t=Sys.time()
full_sim_03 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(3,3),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(0,0)),
  model_fit=list(c('gamma','vm'),c(3,3)),
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

comps[c]=Sys.time()-t
c=c+1

#full_sim_03

################################################################################

# sim1, fit0
t=Sys.time()
full_sim_10 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(0,0),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(1,1)),
  model_fit=list(c('gamma','vm'),c(0,0)),
  N_sim=2,
  N_fit=2,
  n_samples=2000,
  Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),ncol=2),
  delta_sim = c(0.5,0.5),
  param_sim = c(20,40,5,7,
                0,0,2,12),
  autocor_sim = list(matrix(c(0.45,0.55),ncol=1,byrow=TRUE),
                     matrix(c(0.5,0.6),ncol=1,byrow=TRUE)),
  estimate_states = TRUE,
  plot_it = FALSE,
  multicore = TRUE
)

comps[c]=Sys.time()-t
c=c+1

#full_sim_10

################################################################################

# sim1, fit1
t=Sys.time()
full_sim_11 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(1,1),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(1,1)),
  model_fit=list(c('gamma','vm'),c(1,1)),
  N_sim=2,
  N_fit=2,
  n_samples=2000,
  Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),ncol=2),
  delta_sim = c(0.5,0.5),
  param_sim = c(20,40,5,7,
                0,0,2,12),
  autocor_sim = list(matrix(c(0.45,0.55),ncol=1,byrow=TRUE),
                     matrix(c(0.5,0.6),ncol=1,byrow=TRUE)),
  estimate_states = TRUE,
  plot_it = FALSE,
  multicore = TRUE
)

comps[c]=Sys.time()-t
c=c+1

#full_sim_11

################################################################################

# sim1, fit2
t=Sys.time()
full_sim_12 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(2,2),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(1,1)),
  model_fit=list(c('gamma','vm'),c(2,2)),
  N_sim=2,
  N_fit=2,
  n_samples=2000,
  Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),ncol=2),
  delta_sim = c(0.5,0.5),
  param_sim = c(20,40,5,7,
                0,0,2,12),
  autocor_sim = list(matrix(c(0.45,0.55),ncol=1,byrow=TRUE),
                     matrix(c(0.5,0.6),ncol=1,byrow=TRUE)),
  estimate_states = TRUE,
  plot_it = FALSE,
  multicore = TRUE
)

comps[c]=Sys.time()-t
c=c+1

#full_sim_12

################################################################################

# sim1, fit3
t=Sys.time()
full_sim_13 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(3,3),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(1,1)),
  model_fit=list(c('gamma','vm'),c(3,3)),
  N_sim=2,
  N_fit=2,
  n_samples=2000,
  Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),ncol=2),
  delta_sim = c(0.5,0.5),
  param_sim = c(20,40,5,7,
                0,0,2,12),
  autocor_sim = list(matrix(c(0.45,0.55),ncol=1,byrow=TRUE),
                     matrix(c(0.5,0.6),ncol=1,byrow=TRUE)),
  estimate_states = TRUE,
  plot_it = FALSE,
  multicore = TRUE
)

comps[c]=Sys.time()-t
c=c+1

#full_sim_13

################################################################################

# sim2, fit0
t=Sys.time()
full_sim_20 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
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
  plot_it = FALSE,
  multicore = TRUE
)

comps[c]=Sys.time()-t
c=c+1

#full_sim_20

################################################################################

# sim2, fit1
t=Sys.time()
full_sim_21 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(1,1),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(2,2)),
  model_fit=list(c('gamma','vm'),c(1,1)),
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
  plot_it = FALSE,
  multicore = TRUE
)

comps[c]=Sys.time()-t
c=c+1

#full_sim_21

################################################################################

# sim2, fit2
t=Sys.time()
full_sim_22 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
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
  plot_it = FALSE,
  multicore = TRUE
)

comps[c]=Sys.time()-t
c=c+1

#full_sim_22

################################################################################

# sim2, fit3
t=Sys.time()
full_sim_23 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(3,3),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(2,2)),
  model_fit=list(c('gamma','vm'),c(3,3)),
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
  plot_it = FALSE,
  multicore = TRUE
)

comps[c]=Sys.time()-t
c=c+1

#full_sim_23

################################################################################

# sim3, fit0
t=Sys.time()
full_sim_30 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(0,0),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(3,3)),
  model_fit=list(c('gamma','vm'),c(0,0)),
  N_sim=2,
  N_fit=2,
  n_samples=2000,
  Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),ncol=2),
  delta_sim = c(0.5,0.5),
  param_sim = c(20,40,5,7,
                0,0,2,12),
  autocor_sim = list(matrix(c(0.1,0.1,0.25,0.1,0.1,0.35),ncol=3,byrow=TRUE),
                     matrix(c(0.1,0.1,0.3,0.1,0.1,0.4),ncol=3,byrow=TRUE)),
  estimate_states = TRUE,
  plot_it = FALSE,
  multicore = TRUE
)

comps[c]=Sys.time()-t
c=c+1

#full_sim_30

################################################################################

# sim3, fit1
t=Sys.time()
full_sim_31 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(1,1),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(3,3)),
  model_fit=list(c('gamma','vm'),c(1,1)),
  N_sim=2,
  N_fit=2,
  n_samples=2000,
  Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),ncol=2),
  delta_sim = c(0.5,0.5),
  param_sim = c(20,40,5,7,
                0,0,2,12),
  autocor_sim = list(matrix(c(0.1,0.1,0.25,0.1,0.1,0.35),ncol=3,byrow=TRUE),
                     matrix(c(0.1,0.1,0.3,0.1,0.1,0.4),ncol=3,byrow=TRUE)),
  estimate_states = TRUE,
  plot_it = FALSE,
  multicore = TRUE
)

comps[c]=Sys.time()-t
c=c+1

#full_sim_31

################################################################################

# sim3, fit2
t=Sys.time()
full_sim_32 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(2,2),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(3,3)),
  model_fit=list(c('gamma','vm'),c(2,2)),
  N_sim=2,
  N_fit=2,
  n_samples=2000,
  Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),ncol=2),
  delta_sim = c(0.5,0.5),
  param_sim = c(20,40,5,7,
                0,0,2,12),
  autocor_sim = list(matrix(c(0.1,0.1,0.25,0.1,0.1,0.35),ncol=3,byrow=TRUE),
                     matrix(c(0.1,0.1,0.3,0.1,0.1,0.4),ncol=3,byrow=TRUE)),
  estimate_states = TRUE,
  plot_it = FALSE,
  multicore = TRUE
)

comps[c]=Sys.time()-t
c=c+1

#full_sim_32

################################################################################

# sim3, fit3
t=Sys.time()
full_sim_33 <- full_sim_loop(
  simulation = ar_simulation,
  n_runs = 250,
  dists_fitted = c('gamma','vm'),
  p_fitted = c(3,3),
  n_states_fitted = 2,
  n_samples_simulated = 2000,
  # now: simulation parameters
  model_sim=list(c('gamma','vm'),c(3,3)),
  model_fit=list(c('gamma','vm'),c(3,3)),
  N_sim=2,
  N_fit=2,
  n_samples=2000,
  Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),ncol=2),
  delta_sim = c(0.5,0.5),
  param_sim = c(20,40,5,7,
                0,0,2,12),
  autocor_sim = list(matrix(c(0.1,0.1,0.25,0.1,0.1,0.35),ncol=3,byrow=TRUE),
                     matrix(c(0.1,0.1,0.3,0.1,0.1,0.4),ncol=3,byrow=TRUE)),
  estimate_states = TRUE,
  plot_it = FALSE,
  multicore = TRUE
)

comps[c]=Sys.time()-t
c=c+1

#full_sim_33

total_end = Sys.time()
cat("Total simulation time:", total_end-total_begin,"\n")
write.table(comps,
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_updated/computation_times.csv",
            col.names=FALSE, sep=",")


mod=full_sim_33
# find rows where state 1 and 2 are swapped and swap back
swap = which(mod$estimated_parameters$estimated_1_param_1[,1] > 37.5 & mod$estimated_parameters$estimated_1_param_1[,1] < 42.5 & mod$estimated_parameters$estimated_1_param_1[,2] > 17.5 & mod$estimated_parameters$estimated_1_param_1[,2] < 22.5)
for (row in swap){
  hold = mod$estimated_parameters$estimated_1_param_1[row,1]
  mod$estimated_parameters$estimated_1_param_1[row,1]=mod$estimated_parameters$estimated_1_param_1[row,2]
  mod$estimated_parameters$estimated_1_param_1[row,2]=hold
  
  hold = mod$estimated_parameters$estimated_1_param_2[row,1]
  mod$estimated_parameters$estimated_1_param_2[row,1]=mod$estimated_parameters$estimated_1_param_2[row,2]
  mod$estimated_parameters$estimated_1_param_2[row,2]=hold
  
  hold = mod$estimated_parameters$estimated_2_param_1[row,1]
  mod$estimated_parameters$estimated_2_param_1[row,1]=mod$estimated_parameters$estimated_2_param_1[row,2]
  mod$estimated_parameters$estimated_2_param_1[row,2]=hold
  
  hold = mod$estimated_parameters$estimated_2_param_2[row,1]
  mod$estimated_parameters$estimated_2_param_2[row,1]=mod$estimated_parameters$estimated_2_param_2[row,2]
  mod$estimated_parameters$estimated_2_param_2[row,2]=hold
  
  hold = mod$estimated_autocorrelation$estimated_1_autocor[row,1]
  mod$estimated_autocorrelation$estimated_1_autocor[row,1]=mod$estimated_autocorrelation$estimated_1_autocor[row,2]
  mod$estimated_autocorrelation$estimated_1_autocor[row,2]=hold
  
  hold = mod$estimated_autocorrelation$estimated_2_autocor[row,1]
  mod$estimated_autocorrelation$estimated_2_autocor[row,1]=mod$estimated_autocorrelation$estimated_2_autocor[row,2]
  mod$estimated_autocorrelation$estimated_2_autocor[row,2]=hold
  
  mod$decoding_accuracies[row] = 1-mod$decoding_accuracies[row]
}

## save everything

not_global = which(mod$estimated_parameters$estimated_2_param_1[,2] < -0.5 | mod$estimated_parameters$estimated_1_param_1[,2] < 37.5)
not_global

# change filenames according to simulation
write.table(data.frame(c(length(not_global),250-length(not_global)),
                       row.names = c('global optimum not reached','global optimum reached')), 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_updated/full_sim_33/sim_stats.csv",
            col.names=FALSE, sep=",")
write.table(mod$estimated_parameters$estimated_1_param_1, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_updated/full_sim_33/estimated_gamma_mu.csv",
            col.names=FALSE, sep=",")
write.table(mod$estimated_parameters$estimated_1_param_2, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_updated/full_sim_33/estimated_gamma_sigma.csv",
            col.names=FALSE, sep=",")
write.table(mod$estimated_parameters$estimated_2_param_1, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_updated/full_sim_33/estimated_vm_mu.csv",
            col.names=FALSE, sep=",")
write.table(mod$estimated_parameters$estimated_2_param_2, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_updated/full_sim_33/estimated_vm_kappa.csv",
            col.names=FALSE, sep=",")
write.table(mod$estimated_autocorrelation$estimated_1_autocor, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_updated/full_sim_33/estimated_autocor_gamma.csv",
            col.names=FALSE, sep=",")
write.table(mod$estimated_autocorrelation$estimated_2_autocor, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_updated/full_sim_33/estimated_autocor_vm.csv",
            col.names=FALSE, sep=",")
write.table(mod$decoding_accuracies, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_updated/full_sim_33/decoding_accuracies.csv",
            col.names=FALSE, sep=",")











