# 2022-06-05
# Simulations for von Mises HMM

library(MasterThesis)

# data generation
Gamma = matrix(c(0.9,0.1,0.1,0.9),ncol=2,byrow=TRUE)
sample_vm_ar1 <- sample_vonMises_arp(n_samples = 1000,
                    delta=c(0.5,0.5),
                    Gamma=Gamma,
                    N=2,
                    mu=c(-2,2),
                    kappa=c(5,10),
                    autocor=matrix(c(0.4,0.7),nrow=2),
                    p=1)
hist(sample_vm_ar1$data,probability = TRUE)

theta <- c(
  rep(-2,2),
  c(0,0),
  rep(2,2),
  rep(0.1,2)
)
theta.star <- c(
  theta[1:2],
  theta[3:4] * cos(theta[5:6]), # mean
  theta[3:4] * sin(theta[5:6]), # kappa
  qlogis(theta[7:8]) # autocorrelation
)

mllk_vonMises_arp(theta.star, sample_vm_ar1$data, N=2,p=1)
# works

fit_vm_ar0 <- fit_arp_model(mllk_vonMises_arp, sample_vm_ar1$data, theta.star, N=2, 1, 'von Mises')
fit_vm_ar0
# works

## decoding
decoded_states <- viterbi_vonMises_arp(sample_vm_ar1$data, fit_vm_ar0$Gamma, fit_vm_ar0$delta,
                     fit_vm_ar0$autocorrelation,fit_vm_ar0$mu, fit_vm_ar0$kappa,
                     2,0)
sum(decoded_states==sample_vm_ar1$states)/1000



###
### Top to bottom simulation
###

### All in for AR(2) autocorrelation
Gamma_sim = matrix(c(0.8,0.1,0.1,0.1,0.8,0.1,0.1,0.1,0.8),3,3,byrow = TRUE)
Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),2,2,byrow=TRUE)
sim1 <- ar_simulation(model_sim=c('von Mises',1), # autocor simulated model
                      model_fit=c('von Mises',1), # autocor fitted model
                      2, # states simulated model
                      2, # states fitted model
                      1000, # #samples
                      Gamma_sim, # TPM simulated model
                      #delta=c(0.3,0.3,0.4), # Initial distribution simulated model
                      #c(10,20,30), # mu simulated model
                      #c(10,10,10), # sigma simulated model
                      delta_sim=c(0.5,0.5), # Initial distribution simulated model
                      c(0,0,2,10), # mu and kappa simulated model
                      autocor_sim = c(0.1,0.2,0.2,0.3),
                      estimate_states = TRUE,
                      plot_it = TRUE
)
sum(sim1$simulated_model$states==sim1$viterbi_states)/1000
sim1$fitted_model$kappa


##
##
##

#### 
#### now simulations in a loop, 250 for each configuration

n_sims = 250
n_samples = 2000
param_true = c(-0.1,0.1,2,10)
autocor_true = c(0.1,0.1,0.1, 0.2,0.2,0.3)
delta_sim = c(0.5,0.5)
Gamma_true = matrix(c(0.95,0.05,0.05,0.95),2,2,byrow=TRUE)

estimated_mu = matrix(NA,nrow=n_sims,ncol=2)
estimated_kappa = matrix(NA,nrow=n_sims,ncol=2)
estimated_autocor = matrix(NA,nrow=n_sims,ncol=6)
true_states = matrix(NA,nrow=n_sims,ncol=n_samples)
estimated_states = matrix(NA,nrow=n_sims,ncol=n_samples)


## AR(1)
for (i in which(is.na(estimated_mu[,1]))){
  sim <- ar_simulation(model_sim=c('von Mises',3), # autocor simulated model
                       model_fit=c('von Mises',3), # autocor fitted model
                       2, # states simulated model
                       2, # states fitted model
                       n_samples, # #samples
                       Gamma_sim=Gamma_true, # TPM simulated model
                       delta_sim=delta, # Initial distribution simulated model
                       param_sim=param_true, # mu and simulated model
                       autocor_sim = autocor_true,
                       estimate_states = TRUE,
                       plot_it = TRUE
  )
  
  # error handling, skip iteration if optim() in fit function didn't work
  if(anyNA(sim)){
    next
  }
  estimated_mu[i,] = sim$fitted_model$mu
  estimated_kappa[i,] = sim$fitted_model$kappa
  estimated_autocor[i,] = sim$fitted_model$autocorrelation
  true_states[i,] = sim$simulated_model$states
  estimated_states[i,] = sim$viterbi_states
  
}

which(is.na(estimated_mu[,1]))
estimated_mu
estimated_kappa
estimated_autocor
true_states[1:10,1:10]
estimated_states[1:10,1:10]

# find rows where state 1 and 2 are swapped and swap back
swap = which(estimated_mu[,1] > 0 & estimated_mu[,1] < 0.2 & estimated_mu[,2] > -0.2 & estimated_mu[,2] < 0)
for (row in swap){
  hold = estimated_mu[row,1]
  estimated_mu[row,1]=estimated_mu[row,2]
  estimated_mu[row,2]=hold
  
  hold = estimated_kappa[row,1]
  estimated_kappa[row,1]=estimated_kappa[row,2]
  estimated_kappa[row,2]=hold
  
  hold = estimated_autocor[row,1:3]
  estimated_autocor[row,1:3]=estimated_autocor[row,4:6]
  estimated_autocor[row,4:6]=hold
  
  estimated_states[row,]=estimated_states[row,]+1
  estimated_states[row,estimated_states[row,]==3]=1
}

# exclude data where the global optimum is not reached
not_global = which(estimated_kappa[,1] > 4 | estimated_kappa[,2] < 8)
if (length(not_global>0)){ # only delete if there is any to delete
  delete = which(estimated_mu[,1] > 4 | estimated_mu[,2] < 8)
  estimated_mu = estimated_mu[-delete,]
  estimated_kappa = estimated_kappa[-delete,]
  estimated_autocor = estimated_autocor[-delete,]
} 

write.table(data.frame(c(length(not_global),n_sims-length(not_global)),
                       row.names = c('global optimum not reached','global optimum reached')), 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_vonMises_2state_ar3_ar3/sim_stats.csv",
            col.names=FALSE, sep=",")
write.table(estimated_mu, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_vonMises_2state_ar3_ar3/estimated_mu.csv",
            col.names=FALSE, sep=",")
write.table(estimated_kappa, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_vonMises_2state_ar3_ar3/estimated_kappa.csv",
            col.names=FALSE, sep=",")
write.table(estimated_autocor, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_vonMises_2state_ar3_ar3/estimated_autocor.csv",
            col.names=FALSE, sep=",")
write.table(true_states, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_vonMises_2state_ar3_ar3/true_states.csv",
            col.names=FALSE, sep=",")
write.table(estimated_states, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_vonMises_2state_ar3_ar3/estimated_states.csv",
            col.names=FALSE, sep=",")




# [-delete] if necessary
acc = rep(NA,n_sims)#[-delete]
for (i in (1:n_sims)){
  acc[i]=sum(true_states[i,] == estimated_states[i,])/n_samples
}
par(mfrow=c(1,1))
boxplot(acc)



par(mfrow=c(3,2))
boxplot_params(estimated_mu[,1], name=expression(mu[1]), true_value = param_true[1],
               cex.lab=2,cex.axis=1.5)
boxplot_params(estimated_mu[,2], name=expression(mu[2]), true_value = param_true[2],
               cex.lab=2,cex.axis=1.5)
boxplot_params(estimated_kappa[,1], name=expression(kappa[1]), 
               true_value = param_true[3],cex.lab=2,cex.axis=1.5)
boxplot_params(estimated_kappa[,2], name=expression(kappa[2]), 
               true_value = param_true[4],cex.lab=2,cex.axis=1.5)
boxplot_params(estimated_autocor[,1:3], name=expression(phi[1]), 
               true_value = autocor_true[1:3],cex.lab=2,cex.axis=1.5)
boxplot_params(estimated_autocor[,4:6], name=expression(phi[2]), 
               true_value = autocor_true[4:6],cex.lab=2,cex.axis=1.5)

# looking good :)


###
### read in data (if deleted by accident... :D)
###
estimated_mu = read.table(
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_vonMises_2state_ar1_ar1/estimated_mu.csv",
            col.names=FALSE, sep=",")
write.table(estimated_kappa, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_vonMises_2state_ar1_ar1/estimated_kappa.csv",
            col.names=FALSE, sep=",")
write.table(estimated_autocor, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_vonMises_2state_ar1_ar1/estimated_autocor.csv",
            col.names=FALSE, sep=",")
write.table(true_states, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_vonMises_2state_ar1_ar1/true_states.csv",
            col.names=FALSE, sep=",")
write.table(estimated_states, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_vonMises_2state_ar1_ar1/estimated_states.csv",
            col.names=FALSE, sep=",")

