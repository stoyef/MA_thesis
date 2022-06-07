# Master thesis
Thesis code and other things

- Simulation of HMMs with AR(p)-structure in state-dependent process
- Application to real data
- Hopefully much more soon


## ToDo

- [ ] Look over Viterbi function: How can we make sure that the order of the decoded states match the order of the simulated data? 
- [ ] What makes sense for von Mises distribution: Model autocorrelation in $\mu$, $\kappa$ or both?
      ```
      Current code uses only autocorrelation in the parameter mu. For sigma probably another autocorrelation parameter should be used
      ```
- [ ] Alternative plots with densities of AR(0), AR(1), AR(2), AR(3) fit for the same data in one plot?
- [ ] Literature source for the constant constant coefficient of variation (some paper??)
- [x] Re-run of the simulation with 250 runs for updated functions (inclusion of $\sigma$)
- [x] Simulation with 250 runs for von Mises distribution (only autocorrelation in $\mu$), for AR(1), AR(2), AR(3)
- [ ] Incorporate zero-inflation?


## Future prospects

- Implement selected portions of the code (the for loops) in C++ (using the package ```rccp```)
- For real data examples: Incorporate modeling of $\gamma_{ij}$ with covariates? (Not meaningful in simulations, because it distracts from the AR(p) processes)
