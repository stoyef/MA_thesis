# 2022-06-26
# Read in saved simulation statistics and evaluate

sim_stats_33 = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_33/sim_stats.csv",
            row.names = 1,sep=",")
estimated_gamma_mu = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_33/estimated_gamma_mu.csv",
            row.names=1, sep=",")
estimated_gamma_sigma = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_33/estimated_gamma_sigma.csv",
            row.names=1, sep=",")
estimated_vm_mu = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_33/estimated_vm_mu.csv",
            row.names=1, sep=",")
estimated_vm_kappa = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_33/estimated_vm_kappa.csv",
            row.names=1, sep=",")
estimated_autocor_gamma = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_33/estimated_autocor_gamma.csv",
            row.names=1, sep=",")
estimated_autocor_vm = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_33/estimated_autocor_vm.csv",
            row.names=1, sep=",")
decoding_accuracies_20 = read.table("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results/full_sim_20/decoding_accuracies.csv",
            row.names=1, sep=",")


# change according to number of simulation parameters
layout(matrix(c(1,2,3,3),ncol=2,byrow=TRUE))
boxplot(estimated_gamma_mu,xlab=expression(mu),cex.lab=2,xaxt='n',cex.axis=1.5,
        ylim=c(19,43))
axis(1, at=c(1,2), labels=c("",""))
mtext(c(expression(mu[1]),expression(mu[2])),side=1,at=c(1,2),line=1.25,cex=1.5)
abline(h=c(20,40),col=2,lwd=2)
boxplot(estimated_gamma_sigma,xlab=expression(sigma),cex.lab=2,xaxt='n',cex.axis=1.5,
        ylim=c(4,10))
axis(1, at=c(1,2), labels=c("",""))
mtext(c(expression(sigma[1]),expression(sigma[2])),side=1,at=c(1,2),line=1.25,cex=1.5)
abline(h=c(5,7),col=2,lwd=2)
boxplot(estimated_autocor_gamma,xlab=expression(phi),cex.lab=2,cex.axis=1.5,ylim=c(0,1),
        xaxt='n')
#axis(1, at=c(1,2), labels=c("",""))
#mtext(c(expression(phi["1"]),expression(phi["2"])),
#      side=1,at=c(1,2),line=1.25,cex=1.5)
#axis(1, at=c(1,2,3,4), labels=c("","","",""))
#mtext(c(expression(phi["1,1"]),expression(phi["1,2"]),expression(phi["2,1"]),expression(phi["2,2"])),
#      side=1,at=c(1,2,3,4),line=1.25,cex=1.5)
axis(1, at=c(1,2,3,4,5,6), labels=c("","","","","",""))
mtext(c(expression(phi["1,1"]),expression(phi["1,2"]),expression(phi["1,3"]),expression(phi["2,1"]),expression(phi["2,2"]),expression(phi["2,3"])),
      side=1,at=c(1,2,3,4,5,6),line=1.25,cex=1.5)
abline(h=c(0.1,0.25,0.35),col=2,lwd=2) # AR(3)
#abline(h=c(0.15,0.3,0.4),col=2,lwd=2) # AR(2)
#abline(h=c(0.45,0.55),col=2,lwd=2) # AR(1)
#abline(h=0,col=2,lwd=2) # AR(0)
title("Estimated parameters of gamma distribution",outer=TRUE,line=-2,cex.main=2)
#boxplot(estimated_vm_mu,xlab=expression(mu),cex.lab=2,xaxt='n',cex.axis=1.5)
#axis(1, at=c(1,2), labels=c("",""))
#mtext(c(expression(mu[1]),expression(mu[2])),side=1,at=c(1,2),line=1.25,cex=1.5)
#abline(h=c(0),col=2,lwd=2)
#boxplot(estimated_vm_kappa,xlab=expression(kappa),cex.lab=2,xaxt='n',cex.axis=1.5,
#        ylim=c(1,14))
#axis(1, at=c(1,2), labels=c("",""))
#mtext(c(expression(kappa[1]),expression(kappa[2])),side=1,at=c(1,2),line=1.25,cex=1.5)
#abline(h=c(2,12),col=2,lwd=2)
#title("Estimated parameters of von Mises distribution",outer=TRUE,line=-27,cex.main=2)


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
#axis(1, at=c(1,2), labels=c("",""))
#mtext(c(expression(phi["1"]),expression(phi["2"])),
#      side=1,at=c(1,2),line=1.25,cex=1.5)
#axis(1, at=c(1,2,3,4), labels=c("","","",""))
#mtext(c(expression(phi["1,1"]),expression(phi["1,2"]),expression(phi["2,1"]),expression(phi["2,2"])),
#      side=1,at=c(1,2,3,4),line=1.25,cex=1.5)
axis(1, at=c(1,2,3,4,5,6), labels=c("","","","","",""))
mtext(c(expression(phi["1,1"]),expression(phi["1,2"]),expression(phi["1,3"]),expression(phi["2,1"]),expression(phi["2,2"]),expression(phi["2,3"])),
      side=1,at=c(1,2,3,4,5,6),line=1.25,cex=1.5)
#abline(h=0,col=2,lwd=2) # AR(0)
#abline(h=c(0.5,0.6),col=2,lwd=2) # AR(1)
#abline(h=c(0.2,0.3,0.4),col=2,lwd=2) # AR(2)
abline(h=c(0.1,0.3,0.4),col=2,lwd=2) # AR(3)
title("Estimated parameters of von Mises distribution",outer=TRUE,line=-2,cex.main=2)


# Accuracies fit
par(mfrow=c(2,2))
boxplot(c(decoding_accuracies_03, decoding_accuracies_13,
          decoding_accuracies_23, decoding_accuracies_33),ylab="Mean accuracy", 
        main="Accuracies of global decoded states (fitted: AR(3) HMM)",
        bty="n",
        ylim=c(0,1),xaxt='n',xlab='Simulated model')
axis(1, at=c(1,2,3,4), labels=c("","","",""))
mtext(c("AR(0)","AR(1)","AR(2)","AR(3)"),side=1,at=c(1,2,3,4),line=1,cex=1)

# Accuracies sim
#par(mfrow=c(2,2))
boxplot(c(decoding_accuracies_00, decoding_accuracies_01,
          decoding_accuracies_02, decoding_accuracies_03),ylab="Mean accuracy", 
        main="Accuracies of global decoded states (simulated: Normal HMM)",
        bty="n",
        ylim=c(0,1),xaxt='n',xlab='Fitted model')
axis(1, at=c(1,2,3,4), labels=c("","","",""))
mtext(c("AR(0)","AR(1)","AR(2)","AR(3)"),side=1,at=c(1,2,3,4),line=1,cex=1)
boxplot(c(decoding_accuracies_10, decoding_accuracies_11,
          decoding_accuracies_12, decoding_accuracies_13),ylab="Mean accuracy", 
        main="Accuracies of global decoded states (simulated: AR(1) HMM)",
        bty="n",
        ylim=c(0,1),xaxt='n',xlab='Fitted model')
axis(1, at=c(1,2,3,4), labels=c("","","",""))
mtext(c("AR(0)","AR(1)","AR(2)","AR(3)"),side=1,at=c(1,2,3,4),line=1,cex=1)
boxplot(c(decoding_accuracies_20, decoding_accuracies_21,
          decoding_accuracies_22, decoding_accuracies_23),ylab="Mean accuracy", 
        main="Accuracies of global decoded states (simulated: AR(2) HMM)",
        bty="n",
        ylim=c(0,1),xaxt='n',xlab='Fitted model')
axis(1, at=c(1,2,3,4), labels=c("","","",""))
mtext(c("AR(0)","AR(1)","AR(2)","AR(3)"),side=1,at=c(1,2,3,4),line=1,cex=1)
boxplot(c(decoding_accuracies_30, decoding_accuracies_31,
          decoding_accuracies_32, decoding_accuracies_33),ylab="Mean accuracy", 
        main="Accuracies of global decoded states (simulated: AR(3) HMM)",
        bty="n",
        ylim=c(0,1),xaxt='n',xlab='Fitted model')
axis(1, at=c(1,2,3,4), labels=c("","","",""))
mtext(c("AR(0)","AR(1)","AR(2)","AR(3)"),side=1,at=c(1,2,3,4),line=1,cex=1)

for (i in 0:3){
  for (j in 0:3){
    print(paste("Sim:",i, "Fit:", j))
    print(get(paste("sim_stats_",i,j,sep="")))
  }
}
