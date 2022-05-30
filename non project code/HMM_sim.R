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

states_est <- viterbi_arp(x=sim_normal$data, Gamma=mod_normal_data_normal$Gamma, delta=mod_normal_data_normal$delta, 
                          autocor=0, mu=mod_normal_data_normal$mu,
                          sigma=mod_normal_data_normal$sigma, N=2, p=0)
sum(states_est==sim_normal$states)/length(sim_normal$states)

states_est <- viterbi_arp(x=sim_ar1$data, Gamma=mod_ar1_data_ar1$Gamma, delta=mod_ar1_data_ar1$delta, 
                          autocor=mod_ar1_data_ar1$autocorrelation, mu=mod_ar1_data_ar1$mu,
                          sigma=mod_ar1_data_ar1$sigma, N=2, p=1)
sum(states_est==sim_ar1$states)/length(sim_ar1$states)



