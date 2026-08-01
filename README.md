# CATE-SpatioTemporal-Causal

This repository contains the code to reproduce the simulation and application results from the manuscript:

*Estimating Heterogeneous Treatment Effects for Spatio-Temporal Causal Inference*

## Repository Structure
* `Applications/`: Contains scripts and data for real-world application studies.
* `Simulations/`: Includes scripts and configurations for synthetic data simulations.
* `README.md`: This file, providing an overview and instructions.

### Applications Folder

- `Results/`: Saved application results
- `application_code.R`: Script for conducting all empirical analyses
- `app_plot.R`: Script for reproducing application plots

### Simulations Folder

- `Data_specs/`: Data specifications for simulations
- `functions/`: Functions for conducting simulations
- `Iraq_info/`: Iraq window data and distances to main cities/routes
- `Results/`: Saved simulation results
- `SimulationSpatial.R`: Script for running a spatial effect modifier simulation
- `SimulationSpatioTemporal.R`: Script for running a spatio-temporal effect modifier simulation
- `sim_plot.R`: Script for reproducing simulation plots

## Required R Packages

To run the scripts, ensure the following R packages are installed:

- `geocausal` `geocausal` (development version, available at [GitHub repository](https://github.com/lingxiaozhou/geocausal_fork/tree/adaptiveIntervention))
- `spatstat` (Version 3.3-0)
- `tidyverse` (Version 2.0.0)
- `splines` (Version 4.4.0)
- `mgcv` (Version 1.9-1)
- `abind` (Version 1.4-8)

The development version of `geocausal` can be installed directly from GitHub using:

```r
devtools::install_github("lingxiaozhou/geocausal_fork", ref = "adaptiveIntervention")
```




