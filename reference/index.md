# Package index

## Model Fitting

Fit Social Influence Regression models to network data

- [`sir()`](https://netify-dev.github.io/sir/reference/sir.md) : Social
  Influence Regression (SIR) Model
- [`sir_alsfit()`](https://netify-dev.github.io/sir/reference/sir_alsfit.md)
  : Fit SIR Model via Alternating Least Squares (ALS)
- [`sir_optfit()`](https://netify-dev.github.io/sir/reference/sir_optfit.md)
  : Fit SIR Model via Direct Optimization

## S3 Methods

Standard R methods for sir model objects

- [`coef(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/coef.sir.md)
  : Extract Model Coefficients from a SIR Model
- [`vcov(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/vcov.sir.md)
  : Variance-Covariance Matrix for SIR Model Parameters
- [`confint(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/confint.sir.md)
  : Confidence Intervals for SIR Model Parameters
- [`nobs(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/nobs.sir.md)
  : Extract Number of Observations from a SIR Model
- [`fitted(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/fitted.sir.md)
  : Extract Fitted Values from a SIR Model
- [`residuals(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/residuals.sir.md)
  : Extract Residuals from a SIR Model
- [`logLik(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/logLik.sir.md)
  : Extract Log-Likelihood from a SIR Model
- [`AIC(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/AIC.sir.md)
  : Akaike Information Criterion for a SIR Model
- [`BIC(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/BIC.sir.md)
  : Bayesian Information Criterion for a SIR Model
- [`predict(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/predict.sir.md)
  : Predictions from a Fitted SIR Model
- [`summary(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/summary.sir.md)
  : Summary of a Fitted SIR Model
- [`print(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/print.sir.md)
  : Print a Fitted SIR Model
- [`print(`*`<summary.sir>`*`)`](https://netify-dev.github.io/sir/reference/print.summary.sir.md)
  : Print a SIR Model Summary
- [`plot(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/plot.sir.md)
  : Diagnostic Plots for a Fitted SIR Model

## Visualization

Network visualization of influence matrices

- [`plot_sir_network()`](https://netify-dev.github.io/sir/reference/plot_sir_network.md)
  : Network Graph Visualization of Influence Matrices

## Bootstrap & Inference

Bootstrap standard errors and robust inference

- [`boot_sir()`](https://netify-dev.github.io/sir/reference/boot_sir.md)
  : Bootstrap Inference for SIR Model Parameters
- [`print(`*`<boot_sir>`*`)`](https://netify-dev.github.io/sir/reference/print.boot_sir.md)
  : Print Bootstrap SIR Results
- [`summary(`*`<boot_sir>`*`)`](https://netify-dev.github.io/sir/reference/summary.boot_sir.md)
  : Summary of Bootstrap SIR Results
- [`confint(`*`<boot_sir>`*`)`](https://netify-dev.github.io/sir/reference/confint.boot_sir.md)
  : Confidence Intervals from Bootstrap SIR Results

## Counterfactual Scenarios

Build counterfactual scenario arrays for prediction

- [`get_scen_vals()`](https://netify-dev.github.io/sir/reference/get_scen_vals.md)
  : Compute Summary Statistics for Scenario Construction
- [`get_scen_array()`](https://netify-dev.github.io/sir/reference/get_scen_array.md)
  : Build Counterfactual Scenario Array

## Simulation & Data Preparation

Simulate network data and prepare inputs for SIR models

- [`sim_sir()`](https://netify-dev.github.io/sir/reference/sim_sir.md) :
  Simulate Data from a Social Influence Regression Model
- [`rel_covar()`](https://netify-dev.github.io/sir/reference/rel_covar.md)
  : Construct Relational Covariates from a Network Array
- [`cast_array()`](https://netify-dev.github.io/sir/reference/cast_array.md)
  : Cast Directed Dyadic Data into Array Format
- [`eta_tab()`](https://netify-dev.github.io/sir/reference/eta_tab.md) :
  Calculate Linear Predictor (eta) for SIR Model
- [`mll_sir()`](https://netify-dev.github.io/sir/reference/mll_sir.md) :
  Calculate Negative Log-Likelihood for SIR Model

## C++ Backend

Optimized C++ routines for matrix operations

- [`cpp_tprod_A_X_Bt()`](https://netify-dev.github.io/sir/reference/cpp_tprod_A_X_Bt.md)
  : Tensor Product for SIR Model (A \* X \* B')
- [`cpp_amprod_W_v()`](https://netify-dev.github.io/sir/reference/cpp_amprod_W_v.md)
  : Array-Matrix Product for Influence Matrices
- [`cpp_construct_Wbeta_design()`](https://netify-dev.github.io/sir/reference/cpp_construct_Wbeta_design.md)
  : Construct Design Matrix for Alpha Updates in ALS
- [`cpp_construct_Walpha_design()`](https://netify-dev.github.io/sir/reference/cpp_construct_Walpha_design.md)
  : Construct Design Matrix for Beta Updates in ALS
- [`cpp_construct_Wbeta_design_dyn()`](https://netify-dev.github.io/sir/reference/cpp_construct_Wbeta_design_dyn.md)
  : Construct Design Matrix for Alpha Updates with Dynamic W
- [`cpp_construct_Walpha_design_dyn()`](https://netify-dev.github.io/sir/reference/cpp_construct_Walpha_design_dyn.md)
  : Construct Design Matrix for Beta Updates with Dynamic W
- [`cpp_mll_gH()`](https://netify-dev.github.io/sir/reference/cpp_mll_gH.md)
  : Calculate Gradient and Hessian for Direct Optimization
- [`cpp_mll_gH_dyn()`](https://netify-dev.github.io/sir/reference/cpp_mll_gH_dyn.md)
  : Calculate Gradient and Hessian with Dynamic W

## Internal Helpers

Internal utility functions (not exported)

- [`flatten_Y()`](https://netify-dev.github.io/sir/reference/flatten_Y.md)
  : Flatten Y Array for GLM Input
- [`flatten_Z()`](https://netify-dev.github.io/sir/reference/flatten_Z.md)
  : Flatten Z Array for GLM Input
- [`prepare_Z_list()`](https://netify-dev.github.io/sir/reference/prepare_Z_list.md)
  : Prepare Z Array for C++ Consumption
- [`amprod()`](https://netify-dev.github.io/sir/reference/amprod.md) :
  Array-matrix product (R implementation)
- [`mat()`](https://netify-dev.github.io/sir/reference/mat.md) :
  Matricization (R implementation)
- [`tprod()`](https://netify-dev.github.io/sir/reference/tprod.md) :
  Tucker product (R implementation)
