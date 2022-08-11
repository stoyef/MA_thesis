# Master Thesis
Thesis code and other things

- Simulation of HMMs with AR(p)-structure in state-dependent process
- Application to real data
- Hopefully much more soon


## ToDo Theory
- [ ] Literature overview

      - How should I view Markov switching models? Cappe et al (2005) say, Markov switching models are a generalization of HMMs. 
      - This would make my model a Markov switching model and not an HMM...
      - HMMs in general -> Origins of HMMs
      - Markov switching regression 
      - Markov switching autoregression
      - HMMs for animal movement data
      - HMMs for high resolution data (animal movement and other stuff)
- [ ] Theory of HMMs
- [ ] Theory of Markov switching (auto-)regression
- [ ] Model formulation, in general and also in animal movement context
- [ ] Are there any theoretical results for Markov switching regression with AR(p)? Who has done this in the past?
- [ ] Model properties: Number of parameters etc.
- [x] Literature source for the constant coefficient of variation (some paper??) $\to$ Location: Implementation of AR(p) in gamma distribution

## ToDo Simulations
- [x] Look over Viterbi function: How can we make sure that the order of the decoded states match the order of the simulated data? $\to$ Solution: We just check it manually afterwards. Not pretty, but works
- [x] What to do against numerical issues in likelihood computation? If due to bad starting values or unlikely values due to autocorrelation paired with multiplication in multivariate cases there appears a 0 in allprobs, the function return NaN. This leads to a failure in the optimization function... $\to$ This propably is simply a downside. Many errors came from the fact that in theory (due to the weighted mean between autocorrelation and global mean) the $\mu$ value of the gamma distribution could be $<0$. This has been dealt with.
- [x] What makes sense for von Mises distribution: Model autocorrelation in $\mu$, $\kappa$ or both? $\to$ $\mu$ should be modeled in any case, $\kappa$ optionally? $\to$ Modeling $\kappa$ does not make much sense at all. How would we incorporate the autocorrelation? We would have to estimate it from past $\kappa$. This seems kind of weird.
- [x] von Mises distribution: Re-think the weighting of parameter $\mu$ (in the calculation of the Likelihood). Project $[-\pi,\pi]$ onto $\mathbb{R}$ and back? $\to$ Think again if this has already been done in ```starize/unstarize```. Otherwise IMPLEMENT $\to$ Implemented in simulation and density calculation
- [ ] Alternative plots with densities of AR(0), AR(1), AR(2), AR(3) fit for the same data in one plot? IMPLEMENT
- [ ] Comparison of ACFs of simulated and fitted models IMPLEMENT
- [ ] Model selection: Visually, AIC, BIC IMPLEMENT
- [x] Re-run of the gamma distribution simulation with 250 runs for updated functions (inclusion of $\sigma$)
- [x] Simulation with 250 runs for von Mises distribution (only autocorrelation in $\mu$), for AR(1), AR(2), AR(3)
- [x] Master mllk function for arbitrary distribution(s) $\to$ work in progress, how to distinguish parameters that are called the same in multivariate cases???
- [x] Master fit and sample function $\to$ for gamma, von Mises, normal distributions
- [x] Master top to bottom simulation function $\to$ for gamma, von Mises, normal distributions
- [x] Master simulation loop function $\to$ for gamma, von Mises, normal distributions
- [x] Simulation for different configurations: Save results and plots
- [ ] 3D visualization for two-dimensional data IMPLEMENT
- [ ] Think about additional meaningful research prospects for simulations
- [x] Faster computation using multiple cores with ```parallel``` package
- [ ] Re-run simulation with updated data simulation for von Moses distribution

## ToDo Data Example
- [ ] Look for appropriate data

      - One dimensional acceleration data $\to$ gamma distribution
      - Two dimensional step length and turning angle data $\to$ gamma and von Mises distribution
      - Other more sophisticated features
- [ ] Apply different AR(p) HMM models and compare models graphically
- [ ] Model selection: Decide on meaningful criteria/methods


## Future prospects

- Implement selected portions of the code (the for loops) in C++ (using the package ```rccp```)
- Incorporate zero-inflation?
- For real data examples: Incorporate modeling of $\gamma_{ij}$ with covariates? (Not meaningful in simulations, because it distracts from the AR(p) processes)
