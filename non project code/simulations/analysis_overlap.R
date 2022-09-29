# 2022-09-29
# Analysis of the effects of the distribution overlap in step lengths
# ceteris paribus
# A - mu=(20,40)
# B - mu=(25,35)
# C - mu=(25,30)
# A,B,C - sigma=(5,7)


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
    assign(paste("sim_res_2040_",sim,fit,sep=""),
           read_in_sim_data(paste(sim,fit,sep=""), "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_updated/"))
    assign(paste("sim_res_2535_",sim,fit,sep=""),
           read_in_sim_data(paste(sim,fit,sep=""), "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_25_35/"))
    assign(paste("sim_res_2530_",sim,fit,sep=""),
           read_in_sim_data(paste(sim,fit,sep=""), "/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/simulation_results_25_30/"))
  }
}

# Boxplots of estimated step length parameters
# sim: AR(2), fit: AR(0)
par(mfrow=c(1,2))

boxplot(sim_res_2040_20$estimated_gamma_mu[,1], 
        sim_res_2535_20$estimated_gamma_mu[,1], 
        sim_res_2530_20$estimated_gamma_mu[,1],
        pch=19, main=expression(mu[1]),
        xlab="true value"
        )
axis(1, at=c(1,2,3), labels=c("20","25","25"))

boxplot(sim_res_2040_20$estimated_gamma_mu[,2], 
        sim_res_2535_20$estimated_gamma_mu[,2], 
        sim_res_2530_20$estimated_gamma_mu[,2],
        pch=19, main=expression(mu[2]),
        xlab="true value"
        )
axis(1, at=c(1,2,3), labels=c("40","35","30"))


boxplot(sim_res_2040_20$estimated_gamma_sigma[,1], 
        sim_res_2535_20$estimated_gamma_sigma[,1], 
        sim_res_2530_20$estimated_gamma_sigma[,1],
        pch=19, main=expression(sigma[1]),
        xlab="true value"
)
axis(1, at=c(1,2,3), labels=c("5","5","5"))

boxplot(sim_res_2040_20$estimated_gamma_sigma[,2], 
        sim_res_2535_20$estimated_gamma_sigma[,2], 
        sim_res_2530_20$estimated_gamma_sigma[,2],
        pch=19, main=expression(sigma[2]),
        xlab="true value"
)
axis(1, at=c(1,2,3), labels=c("7","7","7"))


# Boxplots of estimated step length parameters
# sim: AR(2), fit: AR(2)
par(mfrow=c(1,2))

boxplot(sim_res_2040_22$estimated_gamma_mu[,1], 
        sim_res_2535_22$estimated_gamma_mu[,1], 
        sim_res_2530_22$estimated_gamma_mu[,1],
        pch=19, main=expression(mu[1]),
        xlab="true value"
)
axis(1, at=c(1,2,3), labels=c("20","25","25"))

boxplot(sim_res_2040_22$estimated_gamma_mu[,2], 
        sim_res_2535_22$estimated_gamma_mu[,2], 
        sim_res_2530_22$estimated_gamma_mu[,2],
        pch=19, main=expression(mu[2]),
        xlab="true value"
)
axis(1, at=c(1,2,3), labels=c("40","35","30"))


boxplot(sim_res_2040_22$estimated_gamma_sigma[,1], 
        sim_res_2535_22$estimated_gamma_sigma[,1], 
        sim_res_2530_22$estimated_gamma_sigma[,1],
        pch=19, main=expression(sigma[1]),
        xlab="true value"
)
axis(1, at=c(1,2,3), labels=c("5","5","5"))

boxplot(sim_res_2040_22$estimated_gamma_sigma[,2], 
        sim_res_2535_22$estimated_gamma_sigma[,2], 
        sim_res_2530_22$estimated_gamma_sigma[,2],
        pch=19, main=expression(sigma[2]),
        xlab="true value"
)
axis(1, at=c(1,2,3), labels=c("7","7","7"))


## Direct comparison 
# Boxplots of estimated step length parameters
# mu
par(mfrow=c(2,4))

boxplot(sim_res_2040_20$estimated_gamma_mu[,1], 
        sim_res_2535_20$estimated_gamma_mu[,1], 
        sim_res_2530_20$estimated_gamma_mu[,1],
        pch=19,# main="simulated: AR(2,2), fitted: AR(0,0)",
        xlab=expression(mu[1]),
        ylim=c(18.75,26.5)
)
axis(1, at=c(1,2,3), labels=c("small","medium","high"),cex.axis=0.8)
points(x=c(1,2,3),y=c(20,25,25),col=2,lwd=1.5,pch=19)

boxplot(sim_res_2040_20$estimated_gamma_mu[,2], 
        sim_res_2535_20$estimated_gamma_mu[,2], 
        sim_res_2530_20$estimated_gamma_mu[,2],
        pch=19, main="",
        xlab=expression(mu[2]),
        ylim=c(27.5,42)
)
axis(1, at=c(1,2,3), labels=c("small","medium","high"),cex.axis=0.8)
points(x=c(1,2,3),y=c(40,35,30),col=2,lwd=1.5,pch=19)
text(-0.5,43.5,"simulated: AR(2,2), fitted: AR(0,0)",xpd='NA',
     font=2,cex=1.3)

boxplot(sim_res_2040_22$estimated_gamma_mu[,1], 
        sim_res_2535_22$estimated_gamma_mu[,1], 
        sim_res_2530_22$estimated_gamma_mu[,1],
        pch=19, #main="simulated: AR(2,2), fitted: AR(2,2)",
        xlab=expression(mu[1]),
        ylim=c(18.75,26.5)
)
axis(1, at=c(1,2,3), labels=c("small","medium","high"),cex.axis=0.8)
points(x=c(1,2,3),y=c(20,25,25),col=2,lwd=1.5,pch=19)

boxplot(sim_res_2040_22$estimated_gamma_mu[,2], 
        sim_res_2535_22$estimated_gamma_mu[,2], 
        sim_res_2530_22$estimated_gamma_mu[,2],
        pch=19, main="",
        xlab=expression(mu[2]),
        ylim=c(27.5,42)
)
text(-0.5,43.5,"simulated: AR(2,2), fitted: AR(2,2)",xpd='NA',
     font=2,cex=1.3)
axis(1, at=c(1,2,3), labels=c("small","medium","high"),cex.axis=0.8)
points(x=c(1,2,3),y=c(40,35,30),col=2,lwd=1.5,pch=19)

# sigma

boxplot(sim_res_2040_20$estimated_gamma_sigma[,1], 
        sim_res_2535_20$estimated_gamma_sigma[,1], 
        sim_res_2530_20$estimated_gamma_sigma[,1],
        pch=19, #main="simulated: AR(2,2), fitted: AR(0,0)",
        xlab=expression(sigma[1]),
        ylim=c(4.25,6.5)
)
axis(1, at=c(1,2,3), labels=c("small","medium","high"),cex.axis=0.8)
points(x=c(1,2,3),y=c(5,5,5),col=2,lwd=1.5,pch=19)

boxplot(sim_res_2040_20$estimated_gamma_sigma[,2], 
        sim_res_2535_20$estimated_gamma_sigma[,2], 
        sim_res_2530_20$estimated_gamma_sigma[,2],
        pch=19, main="",
        xlab=expression(sigma[2]),
        ylim=c(5.75,8.5)
)
axis(1, at=c(1,2,3), labels=c("small","medium","high"),cex.axis=0.8)
points(x=c(1,2,3),y=c(7,7,7),col=2,lwd=1.5,pch=19)


boxplot(sim_res_2040_22$estimated_gamma_sigma[,1], 
        sim_res_2535_22$estimated_gamma_sigma[,1], 
        sim_res_2530_22$estimated_gamma_sigma[,1],
        pch=19, #main="simulated: AR(2,2), fitted: AR(2,2)",
        xlab=expression(sigma[1]),
        ylim=c(4.25,6.5)
)
axis(1, at=c(1,2,3), labels=c("small","medium","high"),cex.axis=0.8)
points(x=c(1,2,3),y=c(5,5,5),col=2,lwd=1.5,pch=19)

boxplot(sim_res_2040_22$estimated_gamma_sigma[,2], 
        sim_res_2535_22$estimated_gamma_sigma[,2], 
        sim_res_2530_22$estimated_gamma_sigma[,2],
        pch=19, main="",
        xlab=expression(sigma[2]),
        ylim=c(5.75,8.5)
)
axis(1, at=c(1,2,3), labels=c("small","medium","high"),cex.axis=0.8)
points(x=c(1,2,3),y=c(7,7,7),col=2,lwd=1.5,pch=19)

legend('topleft',
       c("true values"),
       col=2,
       pch=19,
       lwd=1.5,
       xpd='NA', inset=c(-2.25,-0.25),
       lty=NA)







