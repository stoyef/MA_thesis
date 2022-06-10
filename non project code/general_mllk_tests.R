# test of general mllk function


## von Mises
N=2
p=0
theta <- list(Gamma=rep(-2,2),
              autocor=list(c(0,0)),
              params=list(list(mu=c(-0.5,0.6),
                          kappa=c(8,12))))

theta.star <- list(Gamma=theta$Gamma,
                   autocor=list(qlogis(theta$autocor[[1]])),
                   params=list(list(mu=theta$params[[1]]$mu * cos(theta$params[[1]]$kappa),
                               kappa=theta$params[[1]]$mu * sin(theta$params[[1]]$kappa))))

theta.star

Gamma <- diag(N)
Gamma[!Gamma] <- exp(theta.star$Gamma)
Gamma <- Gamma/rowSums(Gamma)

autocor <- list(exp(theta.star$autocor[[1]]))

# same transformation as in function w2n of moveHMM package:
params <- list(list(mu=Arg(theta.star$params[[1]]$mu+1i*theta.star$params[[1]]$kappa),
               kappa=sqrt(theta.star$params[[1]]$mu^2+theta.star$params[[1]]$kappa^2)))

theta = list(Gamma=Gamma,
             autocor=autocor,
             params=params)


data = sample_vm_ar1$data
mllk(theta=theta, dists='vm', x=data, N=2, p=1)



# 2-dim
theta <- list(Gamma=rep(-2,2),
              autocor=list(c(0,0,0.1,0.1),c(0.2,0.3,0.1,0.4)),
              params=list(list(mu=c(10,20),
                               sigma=c(4,7)),
                          list(mu=c(-1,1),
                               kappa=c(15,12))))

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


