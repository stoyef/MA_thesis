# Master thesis
Thesis code and other things

- Simulation of HMMs with AR(p)-structure in state-dependent process
- Application to real data
- Hopefully much more soon


## ToDo

- [ ] Look over Viterbi function: How can we make sure that the order of the decoded states match the order of the simulated data? 
- [ ] What makes sense for von-Mises distribution: Model autocorrelation in $\mu$, $\kappa$ or both?
      ```html
      // #f03c15
      Current code uses only autocorrelation in the parameter mu. For sigma probably another autocorrelation parameter should be used
      ```
- [ ] Alternative plots with densities of AR(0), AR(1), AR(2), AR(3) fit for the same data in one plot?
- [ ] Source for the constant constant coefficient of variation (some paper??)
- [x] Re-run of the simulation with 250 runs for updated functions (inclusion of sigma)
