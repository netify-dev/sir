# Package index

## Model Fitting

Fit Social Influence Regression models to network data

- [`sir()`](https://netify-dev.github.io/sir/reference/sir.md) : Social
  Influence Regression (SIR) Model

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
- [`plot(`*`<sir_fit>`*`)`](https://netify-dev.github.io/sir/reference/plot.sir_fit.md)
  : Diagnostic Plots for a Fitted SIR Model
- [`tidy(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/tidy.sir.md)
  : Tidy a SIR Model into a Data Frame of Coefficients
- [`glance(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/glance.sir.md)
  : One-Row Model Summary for a SIR Model
- [`augment(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/augment.sir.md)
  : Augment Data With SIR Model Fitted Values and Residuals

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

## Forecasting & Scoring

Forecast held-out networks and score predictions

- [`forecast(`*`<sir_fit>`*`)`](https://netify-dev.github.io/sir/reference/forecast.sir_fit.md)
  [`forecast(`*`<sir>`*`)`](https://netify-dev.github.io/sir/reference/forecast.sir_fit.md)
  : Forecast future networks from a fitted sir model
- [`cv_sir()`](https://netify-dev.github.io/sir/reference/cv_sir.md) :
  Rolling-origin cross-validation for a fitted sir model
- [`print(`*`<sir_cv>`*`)`](https://netify-dev.github.io/sir/reference/print.sir_cv.md)
  : Print rolling-origin cross-validation results
- [`score_sir()`](https://netify-dev.github.io/sir/reference/score_sir.md)
  : Score predicted networks against observed outcomes

## Model-Implied Scenarios

Build scenario arrays for prediction

- [`get_scen_vals()`](https://netify-dev.github.io/sir/reference/get_scen_vals.md)
  : Get Scenario Values for Prediction
- [`get_scen_array()`](https://netify-dev.github.io/sir/reference/get_scen_array.md)
  : Build Scenario Array for Prediction

## Simulation & Data Preparation

Simulate network data and prepare inputs for SIR models

- [`sim_sir()`](https://netify-dev.github.io/sir/reference/sim_sir.md) :
  Simulate Data from a Social Influence Regression Model
- [`icews`](https://netify-dev.github.io/sir/reference/icews.md) : ICEWS
  Inter-State Material Conflict (Monthly)
- [`rel_covar()`](https://netify-dev.github.io/sir/reference/rel_covar.md)
  : Construct Relational Covariates from a Network Array
- [`cast_array()`](https://netify-dev.github.io/sir/reference/cast_array.md)
  : Cast Directed Dyadic Data into Array Format
