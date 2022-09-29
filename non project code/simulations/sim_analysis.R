# 2022-06-26
# Read in saved simulation statistics and evaluate

read_in_sim_data <- function(simfit){ # simfit is a string: 00,01 etc
  sim_stats = read.table(paste("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_25_30/full_sim_",simfit,"/sim_stats.csv",sep=""),
                    row.names = 1,sep=",")
  estimated_gamma_mu = read.table(paste("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_25_30/full_sim_",simfit,"/estimated_gamma_mu.csv",sep=""),
                                  row.names=1, sep=",")
  estimated_gamma_sigma = read.table(paste("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_25_30/full_sim_",simfit,"/estimated_gamma_sigma.csv",sep=""),
                                     row.names=1, sep=",")
  estimated_vm_mu = read.table(paste("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_25_30/full_sim_",simfit,"/estimated_vm_mu.csv",sep=""),
                               row.names=1, sep=",")
  estimated_vm_kappa = read.table(paste("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_25_30/full_sim_",simfit,"/estimated_vm_kappa.csv",sep=""),
                                  row.names=1, sep=",")
  estimated_autocor_gamma = read.table(paste("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_25_30/full_sim_",simfit,"/estimated_autocor_gamma.csv",sep=""),
                                       row.names=1, sep=",")
  estimated_autocor_vm = read.table(paste("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_25_30/full_sim_",simfit,"/estimated_autocor_vm.csv",sep=""),
                                    row.names=1, sep=",")
  decoding_accuracies = read.table(paste("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_25_30/full_sim_",simfit,"/decoding_accuracies.csv",sep=""),
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
    assign(paste("sim_res_",sim,fit,sep=""),
           read_in_sim_data(paste(sim,fit,sep="")))
  }
}

sim_res_03$sim_stats



# fitted: AR(0,0)
boxplot_sim_params(gamma_mu = sim_res_00$estimated_gamma_mu,
                   gamma_sigma = sim_res_00$estimated_gamma_sigma,
                   vm_mu = sim_res_00$estimated_vm_mu,
                   vm_kappa = sim_res_00$estimated_vm_kappa,
                   gamma_auto = 0,
                   vm_auto = 0,
                   true_values = c(20,40,
                                   5,7,
                                   0,0,
                                   2,12),
                   nstates = 2,
                   vm_head_pos=-25)
text = c('true parameters')
legend('topleft',text,col=2,lwd=2,lty=1,#bty='n',
       xpd='NA',inset=c(-0.29,-0.5),
       text.width = 1.05*strwidth(text))

# fitted: AR(1,1)
boxplot_sim_params(gamma_mu = sim_res_31$estimated_gamma_mu,
                   gamma_sigma = sim_res_31$estimated_gamma_sigma,
                   vm_mu = sim_res_31$estimated_vm_mu,
                   vm_kappa = sim_res_31$estimated_vm_kappa,
                   gamma_auto = sim_res_31$estimated_autocor_gamma,
                   vm_auto = sim_res_31$estimated_autocor_vm,
                   true_values = c(20,40,
                                   5,7,
                                   0,0,
                                   2,12,
                                   0.45,0.55,
                                   0.5,0.6),
                   nstates = 2)
text = c('true parameters')
legend('top',text,col=2,lwd=2,lty=1,#bty='n',
       xpd='NA',inset=c(0,-0.4),
       text.width = 1.05*strwidth(text))



# fitted: AR(2,2)
boxplot_sim_params(gamma_mu = sim_res_12$estimated_gamma_mu,
                   gamma_sigma = sim_res_12$estimated_gamma_sigma,
                   vm_mu = sim_res_12$estimated_vm_mu,
                   vm_kappa = sim_res_12$estimated_vm_kappa,
                   gamma_auto = sim_res_12$estimated_autocor_gamma,
                   vm_auto = sim_res_12$estimated_autocor_vm,
                   true_values = c(20,40,
                                   5,7,
                                   0,0,
                                   2,12,
                                   #0,0,0,0,0,0,0,0), #sim:ar(0,0)
                                   0.45,0.55,0.5,0.6), #  sim:ar(1,1)
                                   #0.15,0.3,0.15,0.4,
                                   #0.2,0.3,0.2,0.4),
                   nstates = 2)
text = c('true parameters')
legend('top',text,col=2,lwd=2,lty=1,#bty='n',
       xpd='NA',inset=c(0,-0.4),
       text.width = 1.05*strwidth(text))


# fitted: AR(3,3)
boxplot_sim_params(gamma_mu = sim_res_23$estimated_gamma_mu,
                   gamma_sigma = sim_res_23$estimated_gamma_sigma,
                   vm_mu = sim_res_23$estimated_vm_mu,
                   vm_kappa = sim_res_23$estimated_vm_kappa,
                   gamma_auto = sim_res_23$estimated_autocor_gamma,
                   vm_auto = sim_res_23$estimated_autocor_vm,
                   true_values = c(20,40,
                                   5,7,
                                   0,0,
                                   2,12,
                                   #0,0,0,0,0,0,0,0), #sim:ar(0,0)
                                   #0.45,0.55,0.5,0.6), #  sim:ar(1,1)
                                   #0.15,0.3,0.15,0.4,0.2,0.3,0.2,0.4), # sim:ar(2,2)
                                   0.1,0.1,0.25,0.1,0.1,0.35,
                                   0.1,0.1,0.3,0.1,0.1,0.4),
                   nstates = 2)
text = c('true parameters')
legend('top',text,col=2,lwd=2,lty=1,#bty='n',
       xpd='NA',inset=c(0,-0.4),
       text.width = 1.05*strwidth(text))



# Accuracies fit
par(mfrow=c(2,2))
boxplot(c(sim_res_03$decoding_accuracies, sim_res_13$decoding_accuracies,
          sim_res_23$decoding_accuracies, sim_res_33$decoding_accuracies),ylab="Mean accuracy", 
        main="Accuracies of global decoded states (fitted: AR(3) HMM)",
        bty="n",
        ylim=c(0,1),xaxt='n',xlab='Simulated model')
axis(1, at=c(1,2,3,4), labels=c("","","",""))
mtext(c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),side=1,at=c(1,2,3,4),line=1,cex=1)

# Accuracies sim
par(mfrow=c(2,2))
boxplot(c(sim_res_00$decoding_accuracies, sim_res_01$decoding_accuracies,
          sim_res_02$decoding_accuracies, sim_res_03$decoding_accuracies),ylab="accuracy", 
        main="Simulated: AR(0,0) HMM",
        bty="n",
        ylim=c(0.6,1),xaxt='n',xlab='fitted model',
        pch=19)
axis(1, at=c(1,2,3,4), labels=c("","","",""))
mtext(c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),side=1,at=c(1,2,3,4),line=1,cex=1)
boxplot(c(sim_res_10$decoding_accuracies, sim_res_11$decoding_accuracies,
          sim_res_12$decoding_accuracies, sim_res_13$decoding_accuracies),ylab="accuracy", 
        main="Simulated: AR(1,1) HMM",
        bty="n",
        ylim=c(0.6,1),xaxt='n',xlab='fitted model',
        pch=19)
axis(1, at=c(1,2,3,4), labels=c("","","",""))
mtext(c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),side=1,at=c(1,2,3,4),line=1,cex=1)
boxplot(c(sim_res_20$decoding_accuracies, sim_res_21$decoding_accuracies,
          sim_res_22$decoding_accuracies, sim_res_23$decoding_accuracies),ylab="accuracy", 
        main="Simulated: AR(2,2) HMM",
        bty="n",
        ylim=c(0.6,1),xaxt='n',xlab='fitted model',
        pch=19)
axis(1, at=c(1,2,3,4), labels=c("","","",""))
mtext(c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),side=1,at=c(1,2,3,4),line=1,cex=1)
boxplot(c(sim_res_30$decoding_accuracies, sim_res_31$decoding_accuracies,
          sim_res_32$decoding_accuracies, sim_res_33$decoding_accuracies),ylab="accuracy", 
        main="Simulated: AR(3,3) HMM",
        bty="n",
        ylim=c(0.6,1),xaxt='n',xlab='fitted model',
        pch=19)
axis(1, at=c(1,2,3,4), labels=c("","","",""))
mtext(c("AR(0,0)","AR(1,1)","AR(2,2)","AR(3,3)"),side=1,at=c(1,2,3,4),line=1,cex=1)

for (i in 0:3){
  for (j in 0:3){
    print(paste("Sim:",i, "Fit:", j))
    print(get(paste("sim_res_",i,j,sep=""))$sim_stats)
  }
}
