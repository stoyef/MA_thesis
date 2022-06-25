# 2022-06-26
# Read in saved simulation statistics and evaluate

sim_stats = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_12/sim_stats.csv",
                       row.names = 1,
            col.names=FALSE, sep=",")
estimated_gamma_mu = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_12/estimated_gamma_mu.csv",
            col.names=FALSE, row.names=1, sep=",")
estimated_gamma_sigma = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_12/estimated_gamma_sigma.csv",
            col.names=FALSE, row.names=1, sep=",")
estimated_vm_mu = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_12/estimated_vm_mu.csv",
            col.names=FALSE, row.names=1, sep=",")
estimated_vm_kappa = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_12/estimated_vm_kappa.csv",
            col.names=FALSE, row.names=1, sep=",")
estimated_autocor_gamma = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_12/estimated_autocor_gamma.csv",
            col.names=FALSE, row.names=1, sep=",")
estimated_autocor_vm = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_12/estimated_autocor_vm.csv",
            col.names=FALSE, row.names=1, sep=",")
decoding_accuracies = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_12/decoding_accuracies.csv",
            col.names=FALSE, row.names=1, sep=",")


# change according to number of simulation parameters
layout(matrix(c(1,2,3,3),ncol=2,byrow=TRUE))
boxplot(estimated_gamma_mu,xlab=expression(mu),cex.lab=2,xaxt='n',cex.axis=1.5)
axis(1, at=c(1,2), labels=c("",""))
mtext(c(expression(mu[1]),expression(mu[2])),side=1,at=c(1,2),line=1.25,cex=1.5)
abline(h=c(20,40),col=2,lwd=2)
boxplot(estimated_gamma_sigma,xlab=expression(sigma),cex.lab=2,xaxt='n',cex.axis=1.5)
axis(1, at=c(1,2), labels=c("",""))
mtext(c(expression(sigma[1]),expression(sigma[2])),side=1,at=c(1,2),line=1.25,cex=1.5)
abline(h=c(5,7),col=2,lwd=2)
boxplot(estimated_autocor_gamma,xlab=expression(phi),cex.lab=2,cex.axis=1.5,ylim=c(0,1),
        xaxt='n')
axis(1, at=c(1,2,3,4), labels=c("","","",""))
mtext(c(expression(phi["1,1"]),expression(phi["1,2"]),expression(phi["2,1"]),expression(phi["2,2"])),
      side=1,at=c(1,2,3,4),line=1.25,cex=1.5)
abline(h=autocor_sim[[1]],col=2,lwd=2)
title("Estimated parameters of gamma distribution",outer=TRUE,line=-2,cex.main=2)


layout(matrix(c(1,2,3,3),ncol=2,byrow=TRUE))
boxplot(estimated_vm_mu,xlab=expression(mu),cex.lab=2,xaxt='n',cex.axis=1.5)
axis(1, at=c(1,2), labels=c("",""))
mtext(c(expression(mu[1]),expression(mu[2])),side=1,at=c(1,2),line=1.25,cex=1.5)
abline(h=c(0),col=2,lwd=2)
boxplot(estimated_vm_kappa,xlab=expression(kappa),cex.lab=2,xaxt='n',cex.axis=1.5)
axis(1, at=c(1,2), labels=c("",""))
mtext(c(expression(kappa[1]),expression(kappa[2])),side=1,at=c(1,2),line=1.25,cex=1.5)
abline(h=c(2,12),col=2,lwd=2)
boxplot(estimated_autocor_vm,xlab=expression(phi),cex.lab=2,cex.axis=1.5,ylim=c(0,1),
        xaxt='n')
axis(1, at=c(1,2,3,4), labels=c("","","",""))
mtext(c(expression(phi["1,1"]),expression(phi["1,2"]),expression(phi["2,1"]),expression(phi["2,2"])),
      side=1,at=c(1,2,3,4),line=1.25,cex=1.5)
abline(h=autocor_sim[[2]],col=2,lwd=2)
title("Estimated parameters of von Mises distribution",outer=TRUE,line=-2,cex.main=2)



