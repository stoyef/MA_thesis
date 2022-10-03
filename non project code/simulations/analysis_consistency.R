# 2022-10-03
# Analysis of consistency

## Chapter 3.4.1 Consistency of estimators
## We run the simulation with varying number of samples n to verify consistency
## Choice of n: Can't be much higher than 2000
# n \in \{ 100, 500, 1000, 2000, 5000\}

library(MasterThesis)

# read in data
read_in_sim_data <- function(simfit, path_to_folder){ # simfit is a string: 00,01 etc
  sim_stats = read.table(paste(path_to_folder, "full_sim_",simfit,"/sim_stats.csv",sep=""),
                         row.names = 1,sep=",")
  estimated_gamma_mu = read.table(paste(path_to_folder, "full_sim_",simfit,"/estimated_gamma_mu.csv",sep=""),
                                  row.names=1, sep=",")
  estimated_gamma_sigma = read.table(paste(path_to_folder, "full_sim_",simfit,"/estimated_gamma_sigma.csv",sep=""),
                                     row.names=1, sep=",")
  estimated_vm_mu = read.table(paste(path_to_folder, "full_sim_",simfit,"/estimated_vm_mu.csv",sep=""),
                               row.names=1, sep=",")
  estimated_vm_kappa = read.table(paste(path_to_folder, "full_sim_",simfit,"/estimated_vm_kappa.csv",sep=""),
                                  row.names=1, sep=",")
  estimated_autocor_gamma = read.table(paste(path_to_folder, "full_sim_",simfit,"/estimated_autocor_gamma.csv",sep=""),
                                       row.names=1, sep=",")
  estimated_autocor_vm = read.table(paste(path_to_folder, "full_sim_",simfit,"/estimated_autocor_vm.csv",sep=""),
                                    row.names=1, sep=",")
  decoding_accuracies = read.table(paste(path_to_folder, "full_sim_",simfit,"/decoding_accuracies.csv",sep=""),
                                   row.names=1, sep=",")
  
  ret = list(
    sim_stats = sim_stats,
    estimated_gamma_mu = estimated_gamma_mu,
    estimated_gamma_sigma = estimated_gamma_sigma,
    estimated_vm_mu = estimated_vm_mu,
    estimated_vm_kappa = estimated_vm_kappa,
    estimated_autocor_gamma = estimated_autocor_gamma,
    estimated_autocor_vm = estimated_autocor_vm,
    decoding_accuracies = decoding_accuracies
  )
  return(ret)
}


sim_degs = 0:3
fit_degs = 0:3
for (sim in sim_degs){
  for (fit in fit_degs){
    assign(paste("sim_res_2000_",sim,fit,sep=""),
           read_in_sim_data(paste(sim,fit,sep=""), "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_updated/")) # n=2000
    assign(paste("sim_res_100_",sim,fit,sep=""),
           read_in_sim_data(paste(sim,fit,sep=""), "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_n_100/"))
    assign(paste("sim_res_500_",sim,fit,sep=""),
           read_in_sim_data(paste(sim,fit,sep=""), "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_n_500/"))
    assign(paste("sim_res_1000_",sim,fit,sep=""),
           read_in_sim_data(paste(sim,fit,sep=""), "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_n_1000/"))
    assign(paste("sim_res_5000_",sim,fit,sep=""),
           read_in_sim_data(paste(sim,fit,sep=""), "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_n_5000/"))
  }
}



## consistency of correct models


## step length
# Boxplots of parameters for different n
par(mfrow=c(4,2),
    mai=c(0.5,0.25,0.25,0.25))
boxplot(sim_res_100_22$estimated_gamma_mu[,1],
        sim_res_500_22$estimated_gamma_mu[,1],
        sim_res_2000_22$estimated_gamma_mu[,1],
        sim_res_5000_22$estimated_gamma_mu[,1],
        main=expression(mu[1]),
        pch=19,xlab='',ylim=c(0,53))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=20,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(20,20,20,20),col=2,pch=19)
boxplot(sim_res_100_22$estimated_gamma_mu[,2],
        sim_res_500_22$estimated_gamma_mu[,2],
        sim_res_2000_22$estimated_gamma_mu[,2],
        sim_res_5000_22$estimated_gamma_mu[,2],
        main=expression(mu[2]),
        pch=19,xlab='',ylim=c(0,70))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=40,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(40,40,40,40),col=2,pch=19)
boxplot(sim_res_100_22$estimated_gamma_sigma[,1],
        sim_res_500_22$estimated_gamma_sigma[,1],
        sim_res_2000_22$estimated_gamma_sigma[,1],
        sim_res_5000_22$estimated_gamma_sigma[,1],
        main=expression(sigma[1]),
        pch=19,xlab='', ylim=c(0,14))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=5,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(5,5,5,5),col=2,pch=19)
boxplot(sim_res_100_22$estimated_gamma_sigma[,2],
        sim_res_500_22$estimated_gamma_sigma[,2],
        sim_res_2000_22$estimated_gamma_sigma[,2],
        sim_res_5000_22$estimated_gamma_sigma[,2],
        main=expression(sigma[2]),
        pch=19,xlab='',ylim=c(0,12))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=7,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(7,7,7,7),col=2,pch=19)
legend('topleft', c('true values'), xpd='NA',
       pch=19,col=2, inset=c(-0.16,1.25),
       cex=0.9)

# Boxplots of autoregression params
#par(mfrow=c(2,2),
#    mai=c(0.5,0.5,0.25,0.5))
boxplot(sim_res_100_22$estimated_autocor_gamma[,1],
        sim_res_500_22$estimated_autocor_gamma[,1],
        sim_res_2000_22$estimated_autocor_gamma[,1],
        sim_res_5000_22$estimated_autocor_gamma[,1],
        main=expression(phi["1,1"]),
        pch=19,xlab='',ylim=c(0,1))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=7,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(0.15,0.15,0.15,0.15),col=2,pch=19)
boxplot(sim_res_100_22$estimated_autocor_gamma[,3],
        sim_res_500_22$estimated_autocor_gamma[,3],
        sim_res_2000_22$estimated_autocor_gamma[,3],
        sim_res_5000_22$estimated_autocor_gamma[,3],
        main=expression(phi["2,1"]),
        pch=19,xlab='',ylim=c(0,1))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=7,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(0.15,0.15,0.15,0.15),col=2,pch=19)
boxplot(sim_res_100_22$estimated_autocor_gamma[,2],
        sim_res_500_22$estimated_autocor_gamma[,2],
        sim_res_2000_22$estimated_autocor_gamma[,2],
        sim_res_5000_22$estimated_autocor_gamma[,2],
        main=expression(phi["1,2"]),
        pch=19,xlab='',ylim=c(0,1))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=7,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(0.3,0.3,0.3,0.3),col=2,pch=19)
boxplot(sim_res_100_22$estimated_autocor_gamma[,4],
        sim_res_500_22$estimated_autocor_gamma[,4],
        sim_res_2000_22$estimated_autocor_gamma[,4],
        sim_res_5000_22$estimated_autocor_gamma[,4],
        main=expression(phi["2,2"]),
        pch=19,xlab='',ylim=c(0,1))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=7,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(0.4,0.4,0.4,0.4),col=2,pch=19)
#legend('topleft', c('true values'), xpd='NA',
#       col=2, pch=19, inset=c(-0.22,-0.191),
#       cex=0.8)



## turning angle
# Boxplots of parameters for different n
par(mfrow=c(4,2),
    mai=c(0.5,0.25,0.25,0.25))
boxplot(sim_res_100_22$estimated_vm_mu[,1],
        sim_res_500_22$estimated_vm_mu[,1],
        sim_res_2000_22$estimated_vm_mu[,1],
        sim_res_5000_22$estimated_vm_mu[,1],
        main=expression(mu[1]),
        pch=19,xlab='',ylim=c(-1.5,2.1))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=20,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(0,0,0,0),col=2,pch=19)
boxplot(sim_res_100_22$estimated_vm_mu[,2],
        sim_res_500_22$estimated_vm_mu[,2],
        sim_res_2000_22$estimated_vm_mu[,2],
        sim_res_5000_22$estimated_vm_mu[,2],
        main=expression(mu[2]),
        pch=19,xlab='',ylim=c(-1.5,2.1))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=40,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(0,0,0,0),col=2,pch=19)
boxplot(sim_res_100_22$estimated_vm_kappa[,1],
        sim_res_500_22$estimated_vm_kappa[,1],
        sim_res_2000_22$estimated_vm_kappa[,1],
        sim_res_5000_22$estimated_vm_kappa[,1],
        main=expression(kappa[1]),
        pch=19,xlab='', ylim=c(0,25))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=5,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(2,2,2,2),col=2,pch=19)
boxplot(sim_res_100_22$estimated_vm_kappa[,2],
        sim_res_500_22$estimated_vm_kappa[,2],
        sim_res_2000_22$estimated_vm_kappa[,2],
        sim_res_5000_22$estimated_vm_kappa[,2],
        main=expression(kappa[2]),
        pch=19,xlab='',ylim=c(0,35))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=7,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(12,12,12,12),col=2,pch=19)
legend('topleft', c('true values'), xpd='NA',
       pch=19,col=2, inset=c(-0.16,1.25),
       cex=0.9)

# Boxplots of autoregression params
#par(mfrow=c(2,2),
#    mai=c(0.5,0.5,0.25,0.5))
boxplot(sim_res_100_22$estimated_autocor_vm[,1],
        sim_res_500_22$estimated_autocor_vm[,1],
        sim_res_2000_22$estimated_autocor_vm[,1],
        sim_res_5000_22$estimated_autocor_vm[,1],
        main=expression(phi["1,1"]),
        pch=19,xlab='',ylim=c(0,1))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=7,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(0.2,0.2,0.2,0.2),col=2,pch=19)
boxplot(sim_res_100_22$estimated_autocor_vm[,3],
        sim_res_500_22$estimated_autocor_vm[,3],
        sim_res_2000_22$estimated_autocor_vm[,3],
        sim_res_5000_22$estimated_autocor_vm[,3],
        main=expression(phi["2,1"]),
        pch=19,xlab='',ylim=c(0,1))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=7,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(0.2,0.2,0.2,0.2),col=2,pch=19)
boxplot(sim_res_100_22$estimated_autocor_vm[,2],
        sim_res_500_22$estimated_autocor_vm[,2],
        sim_res_2000_22$estimated_autocor_vm[,2],
        sim_res_5000_22$estimated_autocor_vm[,2],
        main=expression(phi["1,2"]),
        pch=19,xlab='',ylim=c(0,1))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=7,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(0.3,0.3,0.3,0.3),col=2,pch=19)
boxplot(sim_res_100_22$estimated_autocor_vm[,4],
        sim_res_500_22$estimated_autocor_vm[,4],
        sim_res_2000_22$estimated_autocor_vm[,4],
        sim_res_5000_22$estimated_autocor_vm[,4],
        main=expression(phi["2,2"]),
        pch=19,xlab='',ylim=c(0,1))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=7,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(0.4,0.4,0.4,0.4),col=2,pch=19)
#legend('topleft', c('true values'), xpd='NA',
#       col=2, pch=19, inset=c(-0.22,-0.191),
#       cex=0.8)





## step length, AR(0,0) for comparison
# Boxplots of parameters for different n
par(mfrow=c(4,2),
    mai=c(0.5,0.25,0.25,0.25))
boxplot(sim_res_100_00$estimated_gamma_mu[,1],
        sim_res_500_00$estimated_gamma_mu[,1],
        sim_res_2000_00$estimated_gamma_mu[,1],
        sim_res_5000_00$estimated_gamma_mu[,1],
        main=expression(mu[1]),
        pch=19,xlab='',ylim=c(0,53))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=20,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(20,20,20,20),col=2,pch=19)
boxplot(sim_res_100_00$estimated_gamma_mu[,2],
        sim_res_500_00$estimated_gamma_mu[,2],
        sim_res_2000_00$estimated_gamma_mu[,2],
        sim_res_5000_00$estimated_gamma_mu[,2],
        main=expression(mu[2]),
        pch=19,xlab='',ylim=c(0,70))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=40,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(40,40,40,40),col=2,pch=19)
boxplot(sim_res_100_00$estimated_gamma_sigma[,1],
        sim_res_500_00$estimated_gamma_sigma[,1],
        sim_res_2000_00$estimated_gamma_sigma[,1],
        sim_res_5000_00$estimated_gamma_sigma[,1],
        main=expression(sigma[1]),
        pch=19,xlab='', ylim=c(0,14))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=5,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(5,5,5,5),col=2,pch=19)
boxplot(sim_res_100_00$estimated_gamma_sigma[,2],
        sim_res_500_00$estimated_gamma_sigma[,2],
        sim_res_2000_00$estimated_gamma_sigma[,2],
        sim_res_5000_00$estimated_gamma_sigma[,2],
        main=expression(sigma[2]),
        pch=19,xlab='',ylim=c(0,12))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=7,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=c(7,7,7,7),col=2,pch=19)
legend('topleft', c('true values'), xpd='NA',
       pch=19,col=2, inset=c(-0.16,1.25),
       cex=0.9)




## step length
# Boxplots of parameters for different n
# only mu_1 for different ar
par(mfrow=c(2,2),
    mai=c(0.5,0.5,0.25,0.25))
boxplot(sim_res_100_00$estimated_gamma_mu[,2],
        sim_res_500_00$estimated_gamma_mu[,2],
        sim_res_2000_00$estimated_gamma_mu[,2],
        sim_res_5000_00$estimated_gamma_mu[,2],
        main='AR(0,0)-HMM',
        pch=19,xlab='',ylim=c(0,70))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=20,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=rep(40,4),col=2,pch=19)
boxplot(sim_res_100_11$estimated_gamma_mu[,2],
        sim_res_500_11$estimated_gamma_mu[,2],
        sim_res_2000_11$estimated_gamma_mu[,2],
        sim_res_5000_11$estimated_gamma_mu[,2],
        main='AR(1,1)-HMM',
        pch=19,xlab='',ylim=c(0,70))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=40,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=rep(40,4),col=2,pch=19)
boxplot(sim_res_100_22$estimated_gamma_mu[,2],
        sim_res_500_22$estimated_gamma_mu[,2],
        sim_res_2000_22$estimated_gamma_mu[,2],
        sim_res_5000_22$estimated_gamma_mu[,2],
        main='AR(2,2)-HMM',
        pch=19,xlab='', ylim=c(0,70))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=5,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=rep(40,4),col=2,pch=19)
boxplot(sim_res_100_33$estimated_gamma_mu[,2],
        sim_res_500_33$estimated_gamma_mu[,2],
        sim_res_2000_33$estimated_gamma_mu[,2],
        sim_res_5000_33$estimated_gamma_mu[,2],
        main='AR(3,3)-HMM',
        pch=19,xlab='',ylim=c(0,70))
axis(1, at=c(1,2,3,4), labels=c("n = 100","n = 500","n = 2000","n = 5000"),cex.axis=0.9)
#abline(h=7,col=2,lwd=1.5)
points(x=c(1,2,3,4),y=rep(40,4),col=2,pch=19)
legend('topleft', c('true values'), xpd='NA',
       pch=19,col=2, inset=c(-0.23,-0.22),
       cex=0.9)







# Decoding accuracies, maybe for appendix in performance section


for (n in c(100,500,1000,2000,5000)){
  for (i in 0:3){
    for (j in 0:3){
      accs = get(paste('sim_res_',n,'_',i,j,sep=''))$decoding_accuracies
      new_accs = apply(accs,1,function(x) return(max(x,1-x)))
      assign(paste('sim_res_',n,'_',i,j,'$decoding_accuracies',sep=''),
             new_accs)
    }
  }
}

# simulated: AR(2,2), for different n
par(mfrow=c(2,2),
    mai=c(0.5,0.5,0.25,0.5))
boxplot(`sim_res_100_20$decoding_accuracies`,
        `sim_res_100_21$decoding_accuracies`,
        `sim_res_100_22$decoding_accuracies`,
        `sim_res_100_23$decoding_accuracies`,
        pch=19,
        main='n = 100',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_500_20$decoding_accuracies`,
        `sim_res_500_21$decoding_accuracies`,
        `sim_res_500_22$decoding_accuracies`,
        `sim_res_500_23$decoding_accuracies`,
        pch=19,
        main='n = 500',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_1000_20$decoding_accuracies`,
        `sim_res_1000_21$decoding_accuracies`,
        `sim_res_1000_22$decoding_accuracies`,
        `sim_res_1000_23$decoding_accuracies`,
        pch=19,
        main='n = 1000',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_2000_20$decoding_accuracies`,
        `sim_res_2000_21$decoding_accuracies`,
        `sim_res_2000_22$decoding_accuracies`,
        `sim_res_2000_23$decoding_accuracies`,
        pch=19,
        main='n = 2000',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_5000_20$decoding_accuracies`,
        `sim_res_5000_21$decoding_accuracies`,
        `sim_res_5000_22$decoding_accuracies`,
        `sim_res_5000_23$decoding_accuracies`,
        pch=19,
        main='n = 5000',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)



# n=100
par(mfrow=c(2,2))
boxplot(`sim_res_100_00$decoding_accuracies`,
        `sim_res_100_01$decoding_accuracies`,
        `sim_res_100_02$decoding_accuracies`,
        `sim_res_100_03$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(0,0)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_100_10$decoding_accuracies`,
        `sim_res_100_11$decoding_accuracies`,
        `sim_res_100_12$decoding_accuracies`,
        `sim_res_100_13$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(1,1)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_100_20$decoding_accuracies`,
        `sim_res_100_21$decoding_accuracies`,
        `sim_res_100_22$decoding_accuracies`,
        `sim_res_100_23$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(2,2)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_100_30$decoding_accuracies`,
        `sim_res_100_31$decoding_accuracies`,
        `sim_res_100_32$decoding_accuracies`,
        `sim_res_100_33$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(3,3)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)


# n=500
par(mfrow=c(2,2))
boxplot(`sim_res_500_00$decoding_accuracies`,
        `sim_res_500_01$decoding_accuracies`,
        `sim_res_500_02$decoding_accuracies`,
        `sim_res_500_03$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(0,0)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_500_10$decoding_accuracies`,
        `sim_res_500_11$decoding_accuracies`,
        `sim_res_500_12$decoding_accuracies`,
        `sim_res_500_13$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(1,1)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_500_20$decoding_accuracies`,
        `sim_res_500_21$decoding_accuracies`,
        `sim_res_500_22$decoding_accuracies`,
        `sim_res_500_23$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(2,2)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_500_30$decoding_accuracies`,
        `sim_res_500_31$decoding_accuracies`,
        `sim_res_500_32$decoding_accuracies`,
        `sim_res_500_33$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(3,3)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)


# n=1000
par(mfrow=c(2,2))
boxplot(`sim_res_1000_00$decoding_accuracies`,
        `sim_res_1000_01$decoding_accuracies`,
        `sim_res_1000_02$decoding_accuracies`,
        `sim_res_1000_03$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(0,0)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_1000_10$decoding_accuracies`,
        `sim_res_1000_11$decoding_accuracies`,
        `sim_res_1000_12$decoding_accuracies`,
        `sim_res_1000_13$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(1,1)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_1000_20$decoding_accuracies`,
        `sim_res_1000_21$decoding_accuracies`,
        `sim_res_1000_22$decoding_accuracies`,
        `sim_res_1000_23$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(2,2)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_1000_30$decoding_accuracies`,
        `sim_res_1000_31$decoding_accuracies`,
        `sim_res_1000_32$decoding_accuracies`,
        `sim_res_1000_33$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(3,3)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)

# n=2000
par(mfrow=c(2,2))
boxplot(`sim_res_2000_00$decoding_accuracies`,
        `sim_res_2000_01$decoding_accuracies`,
        `sim_res_2000_02$decoding_accuracies`,
        `sim_res_2000_03$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(0,0)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_2000_10$decoding_accuracies`,
        `sim_res_2000_11$decoding_accuracies`,
        `sim_res_2000_12$decoding_accuracies`,
        `sim_res_2000_13$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(1,1)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_2000_20$decoding_accuracies`,
        `sim_res_2000_21$decoding_accuracies`,
        `sim_res_2000_22$decoding_accuracies`,
        `sim_res_2000_23$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(2,2)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_2000_30$decoding_accuracies`,
        `sim_res_2000_31$decoding_accuracies`,
        `sim_res_2000_32$decoding_accuracies`,
        `sim_res_2000_33$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(3,3)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)

# n=5000
par(mfrow=c(2,2))
boxplot(`sim_res_5000_00$decoding_accuracies`,
        `sim_res_5000_01$decoding_accuracies`,
        `sim_res_5000_02$decoding_accuracies`,
        `sim_res_5000_03$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(0,0)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_5000_10$decoding_accuracies`,
        `sim_res_5000_11$decoding_accuracies`,
        `sim_res_5000_12$decoding_accuracies`,
        `sim_res_5000_13$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(1,1)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_5000_20$decoding_accuracies`,
        `sim_res_5000_21$decoding_accuracies`,
        `sim_res_5000_22$decoding_accuracies`,
        `sim_res_5000_23$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(2,2)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)
boxplot(`sim_res_5000_30$decoding_accuracies`,
        `sim_res_5000_31$decoding_accuracies`,
        `sim_res_5000_32$decoding_accuracies`,
        `sim_res_5000_33$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(3,3)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4), labels=c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),cex.axis=0.8)


## Correct vs incorrect model decoding acc
par(mfrow=c(1,2))
boxplot(`sim_res_100_22$decoding_accuracies`,
        `sim_res_500_22$decoding_accuracies`,
        `sim_res_1000_22$decoding_accuracies`,
        `sim_res_2000_22$decoding_accuracies`,
        `sim_res_5000_22$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(2,2)-HMM, fitted: AR(2,2)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4,5), labels=c("n=100","n=500","n=1000","n=2000","n=5000"),cex.axis=0.8)
boxplot(`sim_res_100_20$decoding_accuracies`,
        `sim_res_500_20$decoding_accuracies`,
        `sim_res_1000_20$decoding_accuracies`,
        `sim_res_2000_20$decoding_accuracies`,
        `sim_res_5000_20$decoding_accuracies`,
        pch=19,
        main='Simulated: AR(2,2)-HMM, fitted: AR(0,0)-HMM',
        ylim=c(0.5,1))
axis(1, at=c(1,2,3,4,5), labels=c("n=100","n=500","n=1000","n=2000","n=5000"),cex.axis=0.8)


