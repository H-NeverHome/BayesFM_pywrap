# BayesFM_pywrap
Proudly presenting: a handy python wrapper for the excellent R-package 'BayesFM'. The 'BayesFM'-package is an implementation of the Bayesian version of Exploratory Factor Analysis by [Conti et al. 2014](https://www.sciencedirect.com/science/article/abs/pii/S0304407614001493). It wraps the R-package via rpy2. No knowledge of R required, R-installation and dependency-management handled by wrapper. Going for an sk-learn-ish syntax. 

Functions and methods largely named synonymous in wrapper as  in 'Bayes_FM'-package. For details on function-arguments regarding prior specification, MCMC parameters etc. consult 'Bayes_FM' documentation.

**_This is a very early version. No guarantees of any kind - usage on own risk. Expect bugs. Check results._**

## 'BayesFM' R-Package
- Original Author 'BayesFM': Rémi Piatek 
- ['BayesFM' Documentation & R dependencies](https://cran.r-project.org/web/packages/BayesFM/index.html)
- Implementation of the Bayesian version of Exploratory Factor Analysis by [Conti et al. 2014](https://www.sciencedirect.com/science/article/abs/pii/S0304407614001493)

## Python wrapper 'BayesFM_pywrap' for 'BayesFM'
- Original author: Hauke Niehaus
- Documentatation upcoming
- Minimal use-case upcoming

### Installation
- Wrapper does **not work with Windows**, only unix-systems (will likely not change in the future)  
- Usage of wrapper requires a prior manual R-installation (see [How to Install R on Ubuntu 20.04](https://linuxize.com/post/how-to-install-r-on-ubuntu-20-04/) by Linuxize for help)
- Wrapper not (yet) available as package: Download .zip, extract, move to working directory and import class 'befa' from file into Python-environment

### Minimum Python dependencies for BayesFM_pywrap
- [rpy2](https://rpy2.github.io/)
- [numpy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [arviz (incl. xarray)](https://arviz-devs.github.io/arviz/index.html)

### 'BayesFM' features accessible via BayesFM_pywrap
- Prior simulation/ prior predictive checks
  - Number of active factors
  - Covariance matrix 
-  MCMC fitting procedure with informative & uninformative priors

### 'BayesFM' features (currently) **_not_** accessible via BayesFM_pywrap
- 'simul.dedic.facmod' i.e. generating synthetic data from a dedicated factor model
- Any plotting in R backend (likely adaption to arviz)
- Inclusion of covariates

### (Planned) extensions only available via BayesFM_pywrap
- Plotting via arviz-backend
  - [Autocorrelation-plot](https://arviz-devs.github.io/arviz/api/generated/arviz.plot_autocorr.html#arviz.plot_autocorr)
  - [Trace-plot](https://arviz-devs.github.io/arviz/api/generated/arviz.plot_trace.html)
- Additional MCMC convergence diagnostics (see e.g. ["Convergence Diagnostics for Markov Chain Monte Carlo"](https://www.annualreviews.org/doi/abs/10.1146/annurev-statistics-031219-041300) by Roy(2020))
  - Effective sample size
  - Gelman-Rubin/ R-hat



