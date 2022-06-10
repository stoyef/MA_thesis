# test of general mllk function


## von Mises
N=2
p=0
theta <- list(Gamma=rep(-2,2),
              autocor=list(c(0,0)),
              params=list(list(mu=c(-0.5,0.6),
                          kappa=c(8,12))))

theta.star_vm <- list(Gamma=theta$Gamma,
                   autocor=list(qlogis(theta$autocor[[1]])),
                   params=list(list(mu=theta$params[[1]]$mu * cos(theta$params[[1]]$kappa),
                               kappa=theta$params[[1]]$mu * sin(theta$params[[1]]$kappa))))

theta.star_vm

Gamma <- diag(N)
Gamma[!Gamma] <- exp(theta.star_vm$Gamma)
Gamma <- Gamma/rowSums(Gamma)

autocor <- list(exp(theta.star_vm$autocor[[1]]))

# same transformation as in function w2n of moveHMM package:
params <- list(list(mu=Arg(theta.star_vm$params[[1]]$mu+1i*theta.star_vm$params[[1]]$kappa),
               kappa=sqrt(theta.star_vm$params[[1]]$mu^2+theta.star_vm$params[[1]]$kappa^2)))

theta_vm = list(Gamma=Gamma,
             autocor=autocor,
             params=params)


data_vm = sample_vm_ar1$data
mllk(theta=theta_vm, dists='vm', x=data_vm, N=2, p=1)



# 2-dim
theta <- list(Gamma=rep(-2,2),
              autocor=list(c(0,0,0.1,0.1),c(0.2,0.3,0.1,0.4)),
              params=list(list(mu=c(10,20),
                               sigma=c(4,7)),
                          list(mu=c(-1,1),
                               kappa=c(10,12))))

theta.star <- list(Gamma=theta$Gamma,
                   autocor=list(qlogis(theta$autocor[[1]]),qlogis(theta$autocor[[2]])),
                   params=list(list(mu=log(theta$params[[1]]$mu),
                                    sigma=log(theta$params[[1]]$sigma)),
                               list(mu=theta$params[[2]]$mu * cos(theta$params[[2]]$kappa),
                                    kappa=theta$params[[2]]$mu * sin(theta$params[[2]]$kappa))))

theta.star

Gamma <- diag(N)
Gamma[!Gamma] <- exp(theta.star$Gamma)
Gamma <- Gamma/rowSums(Gamma)

autocor <- list(exp(theta.star$autocor[[1]]), exp(theta.star$autocor[[2]]))

# same transformation as in function w2n of moveHMM package:
params <- list(list(mu=exp(theta.star$params[[1]]$mu),
               sigma=exp(theta.star$params[[1]]$sigma)),
               list(mu=Arg(theta.star$params[[2]]$mu+1i*theta.star$params[[2]]$kappa),
                    kappa=sqrt(theta.star$params[[2]]$mu^2+theta.star$params[[2]]$kappa^2)))

theta = list(Gamma=Gamma,
             autocor=autocor,
             params=params)

data = matrix(c(sim_ar1$data,sample_vm_ar1$data),ncol=2)
mllk(theta=theta, dists=c('gamma', 'vm'), x=data, N=2, p=c(2,2))

# works :)


### now, optimization of mllk


theta.star=c(-2,-2,log(c(10,20,2,4)), 
             c(-1,1)*cos(c(10,12)),c(-1,1)*sin(c(10,12)),qlogis(c(0.1,0.1,0.1,0.1,0.1,0.1)))
mllk(theta.star, dists=c('gamma', 'vm'), x=data, N=2, p=c(1,2))

mod=fit_arp_model(mllk, data, theta.star, N=2, p=c(1,2), dists=c("gamma", "vm"))
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

### 3 dim, 2 state
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


### now, full simulation function

# t.b.c




