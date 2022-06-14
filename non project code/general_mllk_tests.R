# test of general mllk function


## von Mises
N=2
p=1
theta <- c(rep(-2,2),
           c(-0.5,0.6),
           c(8,12),
           c(0,0))

theta.star_vm <- c(theta[1:N*(N-1)],
                   theta[N*(N-1)+1:N] * cos(theta[N*(N-1)+N+1:N]),
                   theta[N*(N-1)+1:N] * sin(theta[N*(N-1)+N+1:N]),
                   qlogis(theta[N*(N-1)+2*N+1:(N*p)]))
theta.star_vm

Gamma_samp = matrix(c(0.9,0.1,0.1,0.9),ncol=2,byrow=TRUE)
sample_vm_ar1 <- sample_vonMises_arp(n_samples = 1000,
                                     delta=c(0.5,0.5),
                                     Gamma=Gamma_samp,
                                     N=2,
                                     mu=c(-2,2),
                                     kappa=c(5,10),
                                     autocor=matrix(c(0.4,0.7),nrow=2),
                                     p=1)
data_vm = sample_vm_ar1$data
mllk(theta.star=theta.star_vm, dists='vm', x=data_vm, N=2, p=1)



# 2-dim
N=2
p=2
theta <- c(rep(-2,2),
           c(10,20,4,5),# Gamma
           c(-1,1,10,12), # von Mises
           c(0.0001,0.0001,0.1,0.1), # autocor gamma
           c(0.2,0.3,0.1,0.4) # autocor von Mises
           )

theta.star_gamma_vm <- c(theta[1:N*(N-1)],
                   log(theta[N*(N-1)+1:(2*N)]), # Gamma
                   theta[N*(N-1)+2*N+1:N] * cos(theta[N*(N-1)+3*N+1:N]), # von Mises
                   theta[N*(N-1)+2*N+1:N] * sin(theta[N*(N-1)+3*N+1:N]), # von Mises
                   qlogis(theta[N*(N-1)+4*N+1:(N*p)]), # autocor gamma
                   qlogis(theta[N*(N-1)+4*N+N*p+1:(N*p)]) # autocor von Mises
                   )
theta.star_gamma_vm

delta_ar1 <- c(0.5,0.5)
Gamma_ar1 <- matrix(c(0.8,0.2,0.2,0.8),nrow=2)
mu_ar1 <- c(10,20)
sigma_ar1 <- c(2,5)
autocor_ar1 <- matrix(c(0.2,0.4),nrow=2)
sim_ar1 <- sample_gamma_arp(1000, delta_ar1, Gamma_ar1, 2, 
                            mu_ar1, sigma_ar1, autocor_ar1, 1)
data = matrix(c(sim_ar1$data,sample_vm_ar1$data),ncol=2)
mllk(theta.star=theta.star_gamma_vm, dists=c('gamma', 'vm'), x=data, N=2, p=c(2,2))

# works :)


### now, optimization of mllk

mllk(theta.star_gamma_vm, dists=c('gamma', 'vm'), x=data, N=2, p=c(2,2))
mod=fit_arp_model(mllk, data, theta.star_gamma_vm, N=2, p=c(2,2), dists=c("gamma", "vm"))
mod
# nice


### now, data sampling in general function

# 2 states, gamma
data = sample_arp(n_samples=1000, 
           delta=c(0.5,0.5),
           Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2),
           N=2,
           params=c(20,40,5,6),
           autocor=0,
           p=0,
           dists=c('gamma'))
data$states
hist(data$data)

# 3 states, gamma
data = sample_arp(n_samples=1000, 
                  delta=c(0.4,0.3,0.3),
                  Gamma=matrix(c(0.9,0.05,0.05,0.05,0.9,0.05,0.05,0.05,0.9),ncol=3,byrow=TRUE),
                  N=3,
                  params=c(20,40,60, 5,6,4),
                  autocor=0,
                  p=0,
                  dists=c('gamma'))
data$states
hist(data$data)

# 2 states, von Mises
data = sample_arp(n_samples=1000, 
                  delta=c(0.5,0.5),
                  Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2),
                  N=2,
                  params=c(-1,1,10,15),
                  autocor=list(matrix(c(0.1,0.2,0.3,0.2),ncol=2,byrow=TRUE)),
                  p=2,
                  dists=c('vm'))
data$states
hist(data$data,breaks=20,prob=T)

# 3 states, von Mises
data = sample_arp(n_samples=1000, 
                  delta=c(0.4,0.3,0.3),
                  Gamma=matrix(c(0.9,0.05,0.05,0.05,0.9,0.05,0.05,0.05,0.9),ncol=3,byrow=TRUE),
                  N=3,
                  params=c(-2,0,2, 10,12,15),
                  autocor=list(matrix(c(0.05,0.05,0.05,0.1,0.2,0.3,0.3,0.2,0.1),ncol=3,byrow=TRUE)),
                  p=3,
                  dists=c('vm'))
data$states
hist(data$data,breaks=30)


### 2 dim, 2 state
data = sample_arp(n_samples=1000, 
                  delta=c(0.5,0.5),
                  Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2),
                  N=2,
                  params=c(-1,1,10,15, -1,1,10,15),
                  autocor=0,#list(matrix(c(0.1,0.2,0.3,0.2),ncol=2,byrow=TRUE)),
                  p=0,
                  dists=c('vm','vm'))
data$states
library(plot3D)
p_x=cut(data$data[,1],30)
p_y=cut(data$data[,2],30)
hist3D(z=table(p_x,p_y))
image2D(z=table(p_x,p_y))

### 2 dim, 2 state
data = sample_arp(n_samples=1000, 
                  delta=c(0.5,0.5),
                  Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2),
                  N=2,
                  params=c(-1,1,10,15, 20,40,2,3.5),
                  autocor=list(matrix(c(0.1,0.2,0.2,0.3),ncol=2,byrow=TRUE),
                               matrix(c(0.2,0.3,0.1,0.1),ncol=2,byrow=TRUE)),
                  p=2,
                  dists=c('vm','gamma'))
data$states
hist(data$data[,1])
hist(data$data[,2])
library(plot3D)
p_x=cut(data$data[,1],30)
p_y=cut(data$data[,2],30)
hist3D(z=table(p_x,p_y))
image2D(z=table(p_x,p_y))


### 3 dim, 2 state
data = sample_arp(n_samples=1000, 
                  delta=c(0.5,0.5),
                  Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2),
                  N=2,
                  params=c(-1,1,10,15, 20,40,2,3.5, -1,2,12,13),
                  autocor=list(matrix(c(0.1,0.2,0.2,0.3),ncol=2,byrow=TRUE),
                               matrix(c(0.2,0.3,0.1,0.1),ncol=2,byrow=TRUE),
                               matrix(c(0.2,0.3,0.1,0.1),ncol=2,byrow=TRUE)),
                  p=2,
                  dists=c('vm','gamma','vm'))
data$states
hist(data$data[,1],breaks=20)
hist(data$data[,2],breaks=20)
hist(data$data[,3],breaks=20)


### now, full simulation function

par(mfrow=c(1,2))
sim_full <- ar_simulation(model_sim=list(c('gamma','vm'),c(1,1)),
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
              estimate_states = FALSE,
              plot_it = TRUE
)


# loop -> simulation without state decoding
n_sims = 250
n_samples = 2000
estimated_gamma_mu = matrix(NA,nrow=n_sims,ncol=2)
estimated_gamma_sigma = matrix(NA,nrow=n_sims,ncol=2)
estimated_vm_mu = matrix(NA,nrow=n_sims,ncol=2)
estimated_vm_kappa = matrix(NA,nrow=n_sims,ncol=2)
estimated_autocor_gamma = matrix(NA,nrow=n_sims,ncol=2)
estimated_autocor_vm = matrix(NA,nrow=n_sims,ncol=2)
#true_states = matrix(NA,nrow=n_sims,ncol=n_samples)
#estimated_states = matrix(NA,nrow=n_sims,ncol=n_samples)

for (i in which(is.na(estimated_gamma_mu[,1]))){
  sim <- ar_simulation(model_sim=list(c('gamma','vm'),c(1,1)),
                                     model_fit=list(c('gamma','vm'),c(1,1)),
                                     N_sim=2,
                                     N_fit=2,
                                     n_samples=n_samples,
                                     Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),ncol=2),
                                     delta_sim = c(0.5,0.5),
                                     param_sim = c(20,40,5,7,
                                                   0,0,2,12),
                                     autocor_sim = list(matrix(c(0.3,0.7),ncol=1),
                                                        matrix(c(0.4,0.7),ncol=1)),
                                     estimate_states = FALSE,
                                     plot_it = TRUE
  )
  
  # error handling, skip iteration if optim() in fit function didn't work
  if(anyNA(sim)){
    next
  }
  estimated_gamma_mu[i,] = sim$fitted_model$params[[1]]$mu
  estimated_gamma_sigma[i,] = sim$fitted_model$params[[1]]$sigma
  estimated_vm_mu[i,] = sim$fitted_model$params[[2]]$mu
  estimated_vm_kappa[i,] = sim$fitted_model$params[[2]]$kappa
  estimated_autocor_gamma[i,] = sim$fitted_model$autocorrelation[[1]]
  estimated_autocor_vm[i,] = sim$fitted_model$autocorrelation[[2]]
  
  #true_states[i,] = sim$simulated_model$states
  #estimated_states[i,] = sim$viterbi_states
  
}

which(is.na(estimated_gamma_mu[,1]))
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






