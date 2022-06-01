## 2022-05-30
## HMM-Simulations

# Read in functions
library(MasterThesis)


#### Goal: Simulate from different kinds of HMMs and compare the performance
####       by refitting 
### - Normal HMM
### - HMM with AR(1)-process in state dependent distribution
### - HMM with higher order AR-processes in state dependent distribution

### We need to specify: Kind of distribution for states, number of states,
###                     Initial distribution, TPM

##
## (Initial) assumption: univariate Gamma distribution, 2 states
##

###
###
### HMM simulation
###
###

### 1. HMM without autocorrelation
# 2 states
delta_normal <- c(0.5,0.5)
Gamma_normal <- matrix(c(0.8,0.2,0.2,0.8),nrow=2)
mu_normal <- c(10,20)
sigma_normal <- c(2,5)
sim_normal <- sample_hmm_normal(1000, delta_normal, Gamma_normal, 2, 
                                mu_normal, sigma_normal)
head(sim_normal$data)

### 2. HMM with AR(1)
# 2 states
delta_ar1 <- c(0.5,0.5)
Gamma_ar1 <- matrix(c(0.8,0.2,0.2,0.8),nrow=2)
mu_ar1 <- c(10,20)
sigma_ar1 <- c(2,5)
autocor_ar1 <- c(0.2,0.4)
sim_ar1 <- sample_hmm_ar1(1000, delta_ar1, Gamma_ar1, 2, 
                          mu_ar1, sigma_ar1, autocor_ar1)
head(sim_ar1$data)

### 3. HMM with AR(2)
# 2 states
delta_ar2 <- c(0.5,0.5)
Gamma_ar2 <- matrix(c(0.8,0.2,0.2,0.8),nrow=2)
mu_ar2 <- c(10,20)
sigma_ar2 <- c(2,5)
autocor_ar2 <- matrix(c(0.05,0.15,0.1,0.3),nrow=2) # first entry for t-2, second entry for t-1
p <- 2
sim_ar2 <- sample_hmm_arp(1000, delta_ar2, Gamma_ar2, 2, 
                          mu_ar2, sigma_ar2, autocor_ar2, p)
head(sim_ar2$data)

### 4. HMM with AR(3)
# 2 states
delta_ar3 <- c(0.5,0.5)
Gamma_ar3 <- matrix(c(0.8,0.2,0.2,0.8),nrow=2)
mu_ar3 <- c(10,20)
sigma_ar3 <- c(2,5)
autocor_ar3 <- matrix(c(0.05,0.05,0.1,0.1,0.1,0.2),nrow=2) # first entry for t-2, second entry for t-1
p <- 3
sim_ar3 <- sample_hmm_arp(1000, delta_ar3, Gamma_ar3, 2, 
                          mu_ar3, sigma_ar3, autocor_ar3, p)
head(sim_ar3$data)


#####
#####
##### Model fitting
#####
#####

### Starting parameters
theta_normal <- c(
  -2,-2, # TPM
  15,30, # mu
  5,10 # sigma
)

theta_normal_star <- c(
  theta_normal[1:2],
  log(theta_normal[3:6])
)

theta_ar1 <- c(
  -2,-2, # TPM
  0.1,0.5, # autocorrelation
  15,30, # mu
  5,10 # sigma
)

theta_ar1_star <- c(
  theta_ar1[1:2],
  qlogis(theta_ar1[3:4]),
  log(theta_ar1[5:8])
)

theta_ar2 <- c(
  -2,-2, # TPM
  0.1,0.1,0.3,0.1, # autocorrelation
  15,30, # mu
  5,10 # sigma
)

theta_ar2_star <- c(
  theta_ar2[1:2],
  qlogis(theta_ar2[3:6]),
  log(theta_ar2[7:10])
)

theta_ar3 <- c(
  -2,-2, # TPM
  0.2,0.05,0.1,0.1,0.05,0.1, # autocorrelation (has to be smaller than 1 in total)
  15,30, # mu
  5,10 # sigma
)

theta_ar3_star <- c(
  theta_ar3[1:2],
  qlogis(theta_ar3[3:8]),
  log(theta_ar3[9:12])
)


##
## Fitting
##

### 1. Data from normal HMM
mod_normal_data_normal <- fit_arp_model(mllk_hmm, sim_normal$data, theta_normal_star, N=2,p=0)
mod_normal_data_normal

mod_ar1_data_normal <- fit_arp_model(mllk_arp, sim_normal$data, theta_ar1_star, N=2,p=1)
mod_ar1_data_normal

mod_ar2_data_normal <- fit_arp_model(mllk_arp, sim_normal$data, theta_ar2_star, N=2,p=2)
mod_ar2_data_normal

mod_ar3_data_normal <- fit_arp_model(mllk_arp, sim_normal$data, theta_ar3_star, N=2,p=3)
mod_ar3_data_normal

### 2. Data from AR(1) HMM
mod_normal_data_ar1 <- fit_arp_model(mllk_hmm, sim_ar1$data, theta_normal_star, N=2,p=0)
mod_normal_data_ar1

mod_ar1_data_ar1 <- fit_arp_model(mllk_arp, sim_ar1$data, theta_ar1_star, N=2,p=1)
mod_ar1_data_ar1

mod_ar2_data_ar1 <- fit_arp_model(mllk_arp, sim_ar1$data, theta_ar2_star, N=2,p=2)
mod_ar2_data_ar1

mod_ar3_data_ar1 <- fit_arp_model(mllk_arp, sim_ar1$data, theta_ar3_star, N=2,p=3)
mod_ar3_data_ar1

### 3. Data from AR(2) HMM
mod_normal_data_ar2 <- fit_arp_model(mllk_hmm, sim_ar2$data, theta_normal_star, N=2,p=0)
mod_normal_data_ar2

mod_ar1_data_ar2 <- fit_arp_model(mllk_arp, sim_ar2$data, theta_ar1_star, N=2,p=1)
mod_ar1_data_ar2

mod_ar2_data_ar2 <- fit_arp_model(mllk_arp, sim_ar2$data, theta_ar2_star, N=2,p=2)
mod_ar2_data_ar2

mod_ar3_data_ar2 <- fit_arp_model(mllk_arp, sim_ar2$data, theta_ar3_star, N=2,p=3)
mod_ar3_data_ar2

### 4. Data from AR(3) HMM
mod_normal_data_ar3 <- fit_arp_model(mllk_hmm, sim_ar3$data, theta_normal_star, N=2,p=0)
mod_normal_data_ar3

mod_ar1_data_ar3 <- fit_arp_model(mllk_arp, sim_ar3$data, theta_ar1_star, N=2,p=1)
mod_ar1_data_ar3

mod_ar2_data_ar3 <- fit_arp_model(mllk_arp, sim_ar3$data, theta_ar2_star, N=2,p=2)
mod_ar2_data_ar3

mod_ar3_data_ar3 <- fit_arp_model(mllk_arp, sim_ar3$data, theta_ar3_star, N=2,p=3)
mod_ar3_data_ar3



###
###
### First model evaluation: AIC and BIC (even though thats not a good idea in general)
###
###

AIC_gamma_HMM(mod_normal_data_normal$mllk,2,0)
AIC_gamma_HMM(mod_ar1_data_normal$mllk,2,1)
AIC_gamma_HMM(mod_ar2_data_normal$mllk,2,2)
AIC_gamma_HMM(mod_ar3_data_normal$mllk,2,3)

AIC_gamma_HMM(mod_normal_data_ar1$mllk,2,0)
AIC_gamma_HMM(mod_ar1_data_ar1$mllk,2,1)
AIC_gamma_HMM(mod_ar2_data_ar1$mllk,2,2)
AIC_gamma_HMM(mod_ar3_data_ar1$mllk,2,3)

AIC_gamma_HMM(mod_normal_data_ar2$mllk,2,0)
AIC_gamma_HMM(mod_ar1_data_ar2$mllk,2,1)
AIC_gamma_HMM(mod_ar2_data_ar2$mllk,2,2)
AIC_gamma_HMM(mod_ar3_data_ar2$mllk,2,3)

AIC_gamma_HMM(mod_normal_data_ar3$mllk,2,0)
AIC_gamma_HMM(mod_ar1_data_ar3$mllk,2,1)
AIC_gamma_HMM(mod_ar2_data_ar3$mllk,2,2)
AIC_gamma_HMM(mod_ar3_data_ar3$mllk,2,3)

# AIC always picks the AR(3)-model

BIC_gamma_HMM(mod_normal_data_normal$mllk,2,0, sim_normal$data)
BIC_gamma_HMM(mod_ar1_data_normal$mllk,2,1, sim_normal$data)
BIC_gamma_HMM(mod_ar2_data_normal$mllk,2,2, sim_normal$data)
BIC_gamma_HMM(mod_ar3_data_normal$mllk,2,3, sim_normal$data)

BIC_gamma_HMM(mod_normal_data_ar1$mllk,2,0, sim_ar1$data)
BIC_gamma_HMM(mod_ar1_data_ar1$mllk,2,1, sim_ar1$data)
BIC_gamma_HMM(mod_ar2_data_ar1$mllk,2,2, sim_ar1$data)
BIC_gamma_HMM(mod_ar3_data_ar1$mllk,2,3, sim_ar1$data)

BIC_gamma_HMM(mod_normal_data_ar2$mllk,2,0, sim_ar2$data)
BIC_gamma_HMM(mod_ar1_data_ar2$mllk,2,1, sim_ar2$data)
BIC_gamma_HMM(mod_ar2_data_ar2$mllk,2,2, sim_ar2$data)
BIC_gamma_HMM(mod_ar3_data_ar2$mllk,2,3, sim_ar2$data)

BIC_gamma_HMM(mod_normal_data_ar3$mllk,2,0, sim_ar3$data)
BIC_gamma_HMM(mod_ar1_data_ar3$mllk,2,1, sim_ar3$data)
BIC_gamma_HMM(mod_ar2_data_ar3$mllk,2,2, sim_ar3$data)
BIC_gamma_HMM(mod_ar3_data_ar3$mllk,2,3, sim_ar3$data)

# BIC is more interesting


###
###
### Global decoding using Viterbi -> for all combination & table in document
###
###

# normal HMM
states_mod_normal_data_normal <- viterbi_arp(x=sim_normal$data, Gamma=mod_normal_data_normal$Gamma, delta=mod_normal_data_normal$delta, 
                          autocor=0, mu=mod_normal_data_normal$mu,
                          sigma=mod_normal_data_normal$sigma, N=2, p=0)
sum(states_mod_normal_data_normal==sim_normal$states)/length(sim_normal$states)

states_mod_normal_data_ar1 <- viterbi_arp(x=sim_ar1$data, Gamma=mod_normal_data_ar1$Gamma, delta=mod_normal_data_ar1$delta, 
                                             autocor=0, mu=mod_normal_data_ar1$mu,
                                             sigma=mod_normal_data_ar1$sigma, N=2, p=0)
sum(states_mod_normal_data_ar1==sim_ar1$states)/length(sim_ar1$states)

states_mod_normal_data_ar2 <- viterbi_arp(x=sim_ar2$data, Gamma=mod_normal_data_ar2$Gamma, delta=mod_normal_data_ar2$delta, 
                                             autocor=0, mu=mod_normal_data_ar2$mu,
                                             sigma=mod_normal_data_ar2$sigma, N=2, p=0)
sum(states_mod_normal_data_ar2==sim_ar2$states)/length(sim_ar2$states)

states_mod_normal_data_ar3 <- viterbi_arp(x=sim_ar3$data, Gamma=mod_normal_data_ar3$Gamma, delta=mod_normal_data_ar3$delta, 
                                          autocor=0, mu=mod_normal_data_ar3$mu,
                                          sigma=mod_normal_data_ar3$sigma, N=2, p=0)
sum(states_mod_normal_data_ar3==sim_ar3$states)/length(sim_ar3$states)

# AR(1)-HMM
states_mod_ar1_data_normal <- viterbi_arp(x=sim_normal$data, Gamma=mod_ar1_data_normal$Gamma, delta=mod_ar1_data_normal$delta, 
                                             autocor=mod_ar1_data_normal$autocorrelation, mu=mod_ar1_data_normal$mu,
                                             sigma=mod_ar1_data_normal$sigma, N=2, p=1)
sum(states_mod_ar1_data_normal==sim_normal$states)/length(sim_normal$states)

states_mod_ar1_data_ar1 <- viterbi_arp(x=sim_ar1$data, Gamma=mod_ar1_data_ar1$Gamma, delta=mod_ar1_data_ar1$delta, 
                                          autocor=mod_ar1_data_ar1$autocorrelation, mu=mod_ar1_data_ar1$mu,
                                          sigma=mod_ar1_data_ar1$sigma, N=2, p=1)
sum(states_mod_ar1_data_ar1==sim_ar1$states)/length(sim_ar1$states)

states_mod_ar1_data_ar2 <- viterbi_arp(x=sim_ar2$data, Gamma=mod_ar1_data_ar2$Gamma, delta=mod_ar1_data_ar2$delta, 
                                          autocor=mod_ar1_data_ar2$autocorrelation, mu=mod_ar1_data_ar2$mu,
                                          sigma=mod_ar1_data_ar2$sigma, N=2, p=1)
sum(states_mod_ar1_data_ar2==sim_ar2$states)/length(sim_ar2$states)

states_mod_ar1_data_ar3 <- viterbi_arp(x=sim_ar3$data, Gamma=mod_ar1_data_ar3$Gamma, delta=mod_ar1_data_ar3$delta, 
                                          autocor=mod_ar1_data_ar3$autocorrelation, mu=mod_ar1_data_ar3$mu,
                                          sigma=mod_ar1_data_ar3$sigma, N=2, p=1)
sum(states_mod_ar1_data_ar3==sim_ar3$states)/length(sim_ar3$states)

# AR(2)-HMM
states_mod_ar2_data_normal <- viterbi_arp(x=sim_normal$data, Gamma=mod_ar2_data_normal$Gamma, delta=mod_ar2_data_normal$delta, 
                                          autocor=mod_ar2_data_normal$autocorrelation, mu=mod_ar2_data_normal$mu,
                                          sigma=mod_ar2_data_normal$sigma, N=2, p=2)
sum(states_mod_ar2_data_normal==sim_normal$states)/length(sim_normal$states)

states_mod_ar2_data_ar1 <- viterbi_arp(x=sim_ar1$data, Gamma=mod_ar2_data_ar1$Gamma, delta=mod_ar2_data_ar1$delta, 
                                       autocor=mod_ar2_data_ar1$autocorrelation, mu=mod_ar2_data_ar1$mu,
                                       sigma=mod_ar2_data_ar1$sigma, N=2, p=2)
sum(states_mod_ar2_data_ar1==sim_ar1$states)/length(sim_ar1$states)

states_mod_ar2_data_ar2 <- viterbi_arp(x=sim_ar2$data, Gamma=mod_ar2_data_ar2$Gamma, delta=mod_ar2_data_ar2$delta, 
                                       autocor=mod_ar2_data_ar2$autocorrelation, mu=mod_ar2_data_ar2$mu,
                                       sigma=mod_ar2_data_ar2$sigma, N=2, p=2)
sum(states_mod_ar2_data_ar2==sim_ar2$states)/length(sim_ar2$states)

states_mod_ar2_data_ar3 <- viterbi_arp(x=sim_ar3$data, Gamma=mod_ar2_data_ar3$Gamma, delta=mod_ar2_data_ar3$delta, 
                                       autocor=mod_ar2_data_ar3$autocorrelation, mu=mod_ar2_data_ar3$mu,
                                       sigma=mod_ar2_data_ar3$sigma, N=2, p=2)
sum(states_mod_ar2_data_ar3==sim_ar3$states)/length(sim_ar3$states)

# AR(3)-HMM
states_mod_ar3_data_normal <- viterbi_arp(x=sim_normal$data, Gamma=mod_ar3_data_normal$Gamma, delta=mod_ar3_data_normal$delta, 
                                          autocor=mod_ar3_data_normal$autocorrelation, mu=mod_ar3_data_normal$mu,
                                          sigma=mod_ar3_data_normal$sigma, N=2, p=3)
sum(states_mod_ar3_data_normal==sim_normal$states)/length(sim_normal$states)

states_mod_ar3_data_ar1 <- viterbi_arp(x=sim_ar1$data, Gamma=mod_ar3_data_ar1$Gamma, delta=mod_ar3_data_ar1$delta, 
                                       autocor=mod_ar3_data_ar1$autocorrelation, mu=mod_ar3_data_ar1$mu,
                                       sigma=mod_ar3_data_ar1$sigma, N=2, p=3)
sum(states_mod_ar3_data_ar1==sim_ar1$states)/length(sim_ar1$states)

states_mod_ar3_data_ar2 <- viterbi_arp(x=sim_ar2$data, Gamma=mod_ar3_data_ar2$Gamma, delta=mod_ar3_data_ar2$delta, 
                                       autocor=mod_ar3_data_ar2$autocorrelation, mu=mod_ar3_data_ar2$mu,
                                       sigma=mod_ar3_data_ar2$sigma, N=2, p=3)
sum(states_mod_ar3_data_ar2==sim_ar2$states)/length(sim_ar2$states)

states_mod_ar3_data_ar3 <- viterbi_arp(x=sim_ar3$data, Gamma=mod_ar3_data_ar3$Gamma, delta=mod_ar3_data_ar3$delta, 
                                       autocor=mod_ar3_data_ar3$autocorrelation, mu=mod_ar3_data_ar3$mu,
                                       sigma=mod_ar3_data_ar3$sigma, N=2, p=3)
sum(states_mod_ar3_data_ar3==sim_ar3$states)/length(sim_ar3$states)


###
### next steps here: take average values of like 10 runs or 1st and 3rd quartiles
### -> better uncertainty quantification
### maybe even confidence intervals???


###
###
### plot the Viterbi decoded states
###
###

plot_states(states_mod_ar3_data_ar3[1:100],names=c("resting","travelling"),
            title=FALSE)

plot_fitted_gamma_dist(data=sim_ar3$data, mu=mod_ar3_data_ar3$mu, sigma=mod_ar3_data_ar3$sigma,
                       delta=mod_ar3_data_ar3$delta,
                       title="2-state-AR(3)-gamma HMM")

plot_data(sim_ar3$data[1:500],name='Step size', title="First 500 observations of AR(3)-gamma HMM")


### All in for no autocorrelation
Gamma_sim = matrix(c(0.8,0.1,0.1,0.1,0.8,0.1,0.1,0.1,0.8),3,3,byrow = TRUE)
Gamma_sim = matrix(c(0.8,0.2,0.2,0.8),2,2,byrow=TRUE)
sim1 <- gamma_simulation(model_sim=0, # autocor simulated model
                 model_fit=0, # autocor fitted model
                 2, # states simulated model
                 2, # states fitted model
                 1000, # #samples
                 Gamma_sim, # TPM simulated model
                 delta=c(0.5,0.5), # Initial distribution simulated model
                 c(10,20), # mu simulated model
                 c(2,5), # sigma simulated model
                 autocor_sim = 0,
                 estimate_states = TRUE,
                 plot_it = TRUE
                 )
sum(sim1$simulated_model$states==sim1$viterbi_states)/1000

### All in for AR(2) autocorrelation
Gamma_sim = matrix(c(0.8,0.1,0.1,0.1,0.8,0.1,0.1,0.1,0.8),3,3,byrow = TRUE)
Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),2,2,byrow=TRUE)
sim1 <- gamma_simulation(model_sim=1, # autocor simulated model
                         model_fit=1, # autocor fitted model
                         2, # states simulated model
                         2, # states fitted model
                         1000, # #samples
                         Gamma_sim, # TPM simulated model
                         #delta=c(0.3,0.3,0.4), # Initial distribution simulated model
                         #c(10,20,30), # mu simulated model
                         #c(10,10,10), # sigma simulated model
                         delta=c(0.5,0.5), # Initial distribution simulated model
                         c(10,30), # mu simulated model
                         c(5,6), # sigma simulated model
                         autocor_sim = c(0.1,0.6),
                         estimate_states = TRUE,
                         plot_it = TRUE
)
sum(sim1$simulated_model$states==sim1$viterbi_states)/1000
sim1$fitted_model$mu


#### 
#### now simulations in a loop, 250 for each configuration

n_sims = 250
n_samples = 2000
mu_true = c(20,40)
sigma_true = c(5,5)
autocor_true = c(0.2,0.6)
delta = c(0.5,0.5)
Gamma_true = matrix(c(0.95,0.05,0.05,0.95),2,2,byrow=TRUE)

estimated_mu = matrix(NA,nrow=n_sims,ncol=2)
estimated_sigma = matrix(NA,nrow=n_sims,ncol=2)
estimated_autocor = matrix(NA,nrow=n_sims,ncol=2)
true_states = matrix(NA,nrow=n_sims,ncol=n_samples)
estimated_states = matrix(NA,nrow=n_sims,ncol=n_samples)


## AR(1)
for (i in which(is.na(estimated_mu[,1]))){
  sim <- gamma_simulation(model_sim=1, # autocor simulated model
                           model_fit=1, # autocor fitted model
                           2, # states simulated model
                           2, # states fitted model
                           n_samples, # #samples
                           Gamma_sim=Gamma_true, # TPM simulated model
                           delta=delta, # Initial distribution simulated model
                           mu_sim=mu_true, # mu simulated model
                           sigma_sim=sigma_true, # sigma simulated model
                           autocor_sim = autocor_true,
                           estimate_states = TRUE,
                           plot_it = TRUE
  )
  
  # error handling, skip iteration if optim() in fit_arp_model() didn't work
  if(anyNA(sim)){
    next
  }
  estimated_mu[i,] = sim$fitted_model$mu
  estimated_sigma[i,] = sim$fitted_model$sigma
  estimated_autocor[i,] = sim$fitted_model$autocorrelation
  true_states[i,] = sim$simulated_model$states
  estimated_states[i,] = sim$viterbi_states
  
}

which(is.na(estimated_mu[,1]))
estimated_mu
estimated_sigma
estimated_autocor
true_states[1:10,1:10]
estimated_states[1:10,1:10]

# exclude data where the global optimum is not reached
not_global = which(estimated_mu[,1] > 25 | estimated_mu[,1] < 15)
if (length(not_global>0)){ # only delete if there is any to delete
  delete = which(estimated_mu[,1] > 25 | estimated_mu[,1] < 15)
  estimated_mu = estimated_mu[-delete,]
  estimated_sigma = estimated_sigma[-delete,]
  estimated_autocor = estimated_autocor[-delete,]
}

write.table(data.frame(c(length(not_global),n_sims-length(not_global)),
                     row.names = c('global optimum not reached','global optimum reached')), 
          "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_100_2state_ar1_ar1/sim_stats.csv",
          col.names=FALSE, sep=",")
write.table(estimated_mu, 
          "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_100_2state_ar1_ar1/estimated_mu.csv",
          col.names=FALSE, sep=",")
write.table(estimated_sigma, 
          "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_100_2state_ar1_ar1/estimated_sigma.csv",
          col.names=FALSE, sep=",")
write.table(estimated_autocor, 
          "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_100_2state_ar1_ar1/estimated_autocor.csv",
          col.names=FALSE, sep=",")
write.table(true_states, 
          "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_100_2state_ar1_ar1/true_states.csv",
          col.names=FALSE, sep=",")
write.table(estimated_states, 
          "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_100_2state_ar1_ar1/estimated_states.csv",
          col.names=FALSE, sep=",")





acc = rep(NA,n_sims)[-delete]
for (i in (1:n_sims)[-delete]){
  acc[i]=sum(true_states[i,] == estimated_states[i,])/n_samples
}
par(mfrow=c(1,1))
boxplot(acc)



par(mfrow=c(3,2))
boxplot_params(estimated_mu[,1], name=expression(mu[1]), true_value = mu_true[1])
boxplot_params(estimated_mu[,2], name=expression(mu[2]), true_value = mu_true[2])
boxplot_params(estimated_sigma[,1], name=expression(sigma[1]), 
               true_value = sigma_true[1])
boxplot_params(estimated_sigma[,2], name=expression(sigma[2]), 
               true_value = sigma_true[2])
boxplot_params(estimated_autocor[,1], name=expression(phi[1]), 
               true_value = autocor_true[1])
boxplot_params(estimated_autocor[,2], name=expression(phi[2]), 
               true_value = autocor_true[2])



