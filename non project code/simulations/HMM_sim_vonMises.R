# 2022-06-05
# Simulations for von Mises HMM

# data generation
Gamma = matrix(c(0.9,0.1,0.1,0.9),ncol=2,byrow=TRUE)
sample_vm_ar1 <- sample_vonMises_arp(n_samples = 1000,
                    delta=c(0.5,0.5),
                    Gamma=Gamma,
                    N=2,
                    mu=c(0,0),
                    kappa=c(5,0.5),
                    autocor=matrix(c(0.2,0.3),nrow=2),
                    p=1)
hist(sample_vm_ar1$data,probability = TRUE)

theta <- c(
  rep(-2,2),
  c(-1,1),
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

fit_vm_ar0 <- fit_arp_model(mllk_vonMises_arp, sample_vm_ar1$data, theta.star, N, 0, 'von Mises')
fit_vm_ar0

curve(fit_vm_ar0$delta[1]*dvm(x,fit_vm_ar0$mu[1],fit_vm_ar0$kappa[1]),-pi,pi,col=2,lwd=1.5,add=TRUE)
curve(fit_vm_ar0$delta[2]*dvm(x,fit_vm_ar0$mu[2],fit_vm_ar0$kappa[2]),-pi,pi,col=4,lwd=1.5,add=TRUE)

