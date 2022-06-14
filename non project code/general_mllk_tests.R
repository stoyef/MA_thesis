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

start_time = Sys.time()
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
  cat(i,'/',n_sims,'\n')
  
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
elapsed_time = Sys.time()-start_time

which(is.na(estimated_gamma_mu[,1]))
estimated_gamma_mu
estimated_gamma_sigma
estimated_vm_mu
estimated_vm_kappa
estimated_autocor_gamma
estimated_autocor_vm

# find rows where state 1 and 2 are swapped and swap back
swap = which(estimated_gamma_mu[,1] > 35 & estimated_gamma_mu[,1] < 45 & estimated_gamma_mu[,2] > 15 & estimated_gamma_mu[,2] < 25)
for (row in swap){
  hold = estimated_gamma_mu[row,1]
  estimated_gamma_mu[row,1]=estimated_gamma_mu[row,2]
  estimated_gamma_mu[row,2]=hold
  
  hold = estimated_gamma_sigma[row,1]
  estimated_gamma_sigma[row,1]=estimated_gamma_sigma[row,2]
  estimated_gamma_sigma[row,2]=hold
  
  hold = estimated_vm_mu[row,1]
  estimated_vm_mu[row,1]=estimated_vm_mu[row,2]
  estimated_vm_mu[row,2]=hold
  
  hold = estimated_vm_kappa[row,1]
  estimated_vm_kappa[row,1]=estimated_vm_kappa[row,2]
  estimated_vm_kappa[row,2]=hold
  
  hold = estimated_autocor_gamma[row,1]
  estimated_autocor_gamma[row,1]=estimated_autocor_gamma[row,2]
  estimated_autocor_gamma[row,2]=hold
  
  hold = estimated_autocor_vm[row,1]
  estimated_autocor_vm[row,1]=estimated_autocor_vm[row,2]
  estimated_autocor_vm[row,2]=hold
  
#  estimated_states[row,]=estimated_states[row,]+1
#  estimated_states[row,estimated_states[row,]==3]=1
}

# exclude data where the global optimum is not reached
not_global = which(estimated_gamma_mu[,1] > 25 | estimated_gamma_mu[,2] < 35)
if (length(not_global>0)){ # only delete if there is any to delete
  delete = which(estimated_gamma_mu[,1] > 25 | estimated_gamma_mu[,2] < 35)
  estimated_gamma_mu = estimated_gamma_mu[-delete,]
  estimated_gamma_sigma = estimated_gamma_sigma[-delete,]
  estimated_vm_mu = estimated_vm_mu[-delete,]
  estimated_vm_kappa = estimated_vm_kappa[-delete,]
  estimated_autocor_gamma = estimated_autocor_gamma[-delete,]
  estimated_autocor_vm = estimated_autocor_vm[-delete,]
  } 

write.table(data.frame(c(length(not_global),n_sims-length(not_global)),
                       row.names = c('global optimum not reached','global optimum reached')), 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_gamma_vm_2state_ar1_ar1/sim_stats.csv",
            col.names=FALSE, sep=",")
write.table(estimated_gamma_mu, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_gamma_vm_2state_ar1_ar1/estimated_gamma_mu.csv",
            col.names=FALSE, sep=",")
write.table(estimated_gamma_sigma, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_gamma_vm_2state_ar1_ar1/estimated_gamma_sigma.csv",
            col.names=FALSE, sep=",")
write.table(estimated_vm_mu, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_gamma_vm_2state_ar1_ar1/estimated_vm_mu.csv",
            col.names=FALSE, sep=",")
write.table(estimated_vm_kappa, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_gamma_vm_2state_ar1_ar1/estimated_vm_kappa.csv",
            col.names=FALSE, sep=",")
write.table(estimated_autocor_gamma, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_gamma_vm_2state_ar1_ar1/estimated_autocor_gamma.csv",
            col.names=FALSE, sep=",")
write.table(estimated_autocor_vm, 
            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_gamma_vm_2state_ar1_ar1/estimated_autocor_vm.csv",
            col.names=FALSE, sep=",")
#write.table(true_states, 
#            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_gamma_vm_2state_ar1_ar1/true_states.csv",
#            col.names=FALSE, sep=",")
#write.table(estimated_states, 
#            "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation results/sim_250_gamma_vm_2state_ar1_ar1/estimated_states.csv",
#            col.names=FALSE, sep=",")


# [-delete] if necessary
#acc = rep(NA,n_sims)#[-delete]
#for (i in (1:n_sims)){
#  acc[i]=sum(true_states[i,] == estimated_states[i,])/n_samples
#}
#par(mfrow=c(1,1))
#boxplot(acc)



par(mfrow=c(3,2))
boxplot_params(estimated_gamma_mu[,1], name=expression(mu[1]), true_value = param_sim[1])
boxplot_params(estimated_gamma_mu[,2], name=expression(mu[2]), true_value = param_sim[2])
boxplot_params(estimated_gamma_sigma[,1], name=expression(sigma[1]), 
               true_value = param_sim[3])
boxplot_params(estimated_gamma_sigma[,2], name=expression(sigma[2]), 
               true_value = param_sim[4])
boxplot_params(estimated_autocor_gamma[,1], name=expression(phi[1]), 
               true_value = autocor_sim[[1]][1])
boxplot_params(estimated_autocor_gamma[,2], name=expression(phi[2]), 
               true_value = autocor_sim[[1]][2])

par(mfrow=c(1,2))
boxplot(estimated_gamma_mu,ylim=c(19,42),xlab=expression(mu))
abline(h=c(20,40),col=2,lwd=1.5)
boxplot(estimated_gamma_sigma,xlab=expression(sigma))
abline(h=c(5,7),col=2,lwd=1.5)
title("Parameters of the gamma distribution",outer=TRUE,line=-3)

boxplot(estimated_vm_mu,xlab=expression(mu))
abline(h=c(0),col=2,lwd=1.5)
boxplot(estimated_vm_kappa,xlab=expression(kappa))
abline(h=c(2,12),col=2,lwd=1.5)
title("Parameters of the von Mises distribution",outer=TRUE,line=-3)
