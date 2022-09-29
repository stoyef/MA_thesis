# Master Thesis
Thesis code (structured as an R-package). Additional code to generate simulations, data applications and additional tests are located in [non project code]

Functionalities include:
- Computation of AR(p)-HMM likelihood
- Fitting AR(p)-HMMs
- Global decoding
- Pseudo residuals
- Simulation of data from AR(p)-HMMs
- Applications to real data


## ToDo Theory
- [x] Literature overview
- [x] Theory of HMMs
- [x] Theory of Markov switching (auto-)regression
- [x] Model formulation, in general and also in animal movement context
- [x] Are there any theoretical results for Markov switching regression with AR(p)? Who has done this in the past?
- [x] Model properties: Number of parameters etc.
- [x] Literature source for the constant coefficient of variation (some paper??) $\to$ Location: Implementation of AR(p) in gamma distribution
- [ ] Better ending for chapter 2, maybe a short summary?

## ToDo Simulations
- [x] Look over Viterbi function: How can we make sure that the order of the decoded states match the order of the simulated data? $\to$ Solution: We just check it manually afterwards. Not pretty, but works
- [x] What to do against numerical issues in likelihood computation? If due to bad starting values or unlikely values due to autocorrelation paired with multiplication in multivariate cases there appears a 0 in allprobs, the function return NaN. This leads to a failure in the optimization function... $\to$ This propably is simply a downside. Many errors came from the fact that in theory (due to the weighted mean between autocorrelation and global mean) the $\mu$ value of the gamma distribution could be $<0$. This has been dealt with. CORRECTION: The errors were generated by a bug (a not well thought through calculation of starting values for the autoregressive parameters could lead to a sum of >1 in those parameters which could lead to a value <0 for the mean value)
- [x] What makes sense for von Mises distribution: Model autocorrelation in $\mu$, $\kappa$ or both? $\to$ $\mu$ should be modeled in any case, $\kappa$ optionally? $\to$ Modeling $\kappa$ does not make much sense at all. How would we incorporate the autocorrelation? We would have to estimate it from past $\kappa$. This seems kind of weird
- [x] von Mises distribution: Re-think the weighting of parameter $\mu$ (in the calculation of the Likelihood). Project $[-\pi,\pi]$ onto $\mathbb{R}$ and back? $\to$ Think again if this has already been done in ```starize/unstarize```. Otherwise IMPLEMENT $\to$ Implemented in simulation and density calculation
- [ ] Alternative plots with densities of AR(0), AR(1), AR(2), AR(3) fit for the same data in one plot? IMPLEMENT
- [ ] Model selection: AIC, BIC IMPLEMENT
- [ ] Model evaluation using Pseudo residuals? $\to$ Only if this adds interesting info, otherwise too much
- [ ] Re-run of everything, before writing :)
- [x] Re-run of the gamma distribution simulation with 250 runs for updated functions (inclusion of $\sigma$)
- [x] Simulation with 250 runs for von Mises distribution (only autocorrelation in $\mu$), for AR(1), AR(2), AR(3)
- [x] Master mllk function for arbitrary distribution(s) $\to$ work in progress, how to distinguish parameters that are called the same in multivariate cases???
- [x] Master fit and sample function $\to$ for gamma, von Mises, normal distributions
- [x] Master top to bottom simulation function $\to$ for gamma, von Mises, normal distributions
- [x] Master simulation loop function $\to$ for gamma, von Mises, normal distributions
- [x] Simulation for different configurations: Save results and plots
- [ ] Think about additional meaningful research prospects for simulations
- [x] Faster computation using multiple cores with ```parallel``` package
- [ ] Table for computation speed should display relative speed compared to normal HMM to avoid user bias (absolute values in appendix)

## ToDo Data Example
- [x] Look for appropriate data
- [x] Apply different AR(p) HMM models and compare models graphically
- [ ] Parametric Bootstrap for confidence intervals of the parameters? Would be extremely time consuming. Maybe only for one example? (Technique Re-fit Bootstrap samples that are generated from HMM with MLE parameters, do this a lot of times)
- [x] Model selection: Decide on meaningful criteria/methods
- [x] Check models with pseudo residuals $\to$ Progress in goodness of fit can be displayed here $\to$ look up how to compute them in source code of moveHMM

## ToDo At the end
- [ ] Look through thesis with focus on dimension. Is the bold face for vectors consistently used? Write at the first time it happens that we most often only handle the univariate case
- [ ] Hut über alle geschätzten Parameter
- [ ] autoregressIVE parameters

## Future prospects

- Implement selected portions of the code (the for loops) in C++ (using the package ```rccp```)
- Incorporate zero-inflation?
- For real data examples: Incorporate modeling of $\gamma_{ij}$ with covariates? (Not meaningful in simulations, because it distracts from the AR(p) processes)
