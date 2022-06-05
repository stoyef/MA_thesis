# 2022-06-03
# More simulations of Gamma HMMs with autocorrelation

# Read in functions
library(MasterThesis)

#### Include sigma in modeling of autocorrelation, using constant coefficient 
#### of variation

### All in for AR(1) autocorrelation
Gamma_sim = matrix(c(0.8,0.1,0.1,0.05,0.9,0.05,0.01,0.04,0.95),3,3,byrow = TRUE)
#Gamma_sim = matrix(c(0.9,0.1,0.1,0.9),2,2,byrow=TRUE)
sim1 <- gamma_simulation(model_sim=1, # autocor simulated model
                            model_fit=1, # autocor fitted model
                            3, # states simulated model
                            3, # states fitted model
                            2500, # #samples
                            Gamma_sim, # TPM simulated model
                            delta=c(0.4,0.3,0.3), # Initial distribution simulated model
                            c(10,20,30), # mu simulated model
                            c(5,6,4), # sigma simulated model
                            #delta=c(0.5,0.5), # Initial distribution simulated model
                            #c(10,30), # mu simulated model
                            #c(5,6), # sigma simulated model
                            autocor_sim = c(0.1,0.6,0.9),
                            estimate_states = TRUE,
                            plot_it = TRUE
)
sum(sim1$simulated_model$states==sim1$viterbi_states)/2500
sim1$fitted_model$mu
sim1$fitted_model$delta
