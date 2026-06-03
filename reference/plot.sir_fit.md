# Diagnostic Plots for a Fitted SIR Model

Produces a selection of diagnostic plots for model assessment. By
default, plots 1-4 (influence matrix heatmaps and distributions) are
shown. Use the `which` argument to select specific plots. All plots use
`ggplot2` and are combined via `patchwork` when `combine = TRUE`.

## Usage

``` r
# S3 method for class 'sir_fit'
plot(
  x,
  which = 1:4,
  combine = TRUE,
  title = NULL,
  node_labels = NULL,
  period = NULL,
  theme_base = theme_bw(),
  ...
)
```

## Arguments

- x:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- which:

  Integer vector selecting which plots to produce. Options:

  1

  :   Heatmap of sender influence matrix A. Shows how each node
      influences others' outgoing ties.

  2

  :   Heatmap of receiver influence matrix B. Shows how each node
      affects others' incoming ties.

  3

  :   Histogram and density of off-diagonal A values. Useful for
      assessing the overall strength and distribution of sender effects.

  4

  :   Histogram and density of off-diagonal B values. Same for receiver
      effects.

  5

  :   Convergence trace plot showing deviance across ALS iterations.
      Iteration history is always stored in the fitted model.

  6

  :   Coefficient plot with 95% confidence intervals. Requires standard
      errors (`calc_se = TRUE`). Parameters are grouped by type
      (exogenous, alpha, beta).

- combine:

  Logical. If TRUE (default), combines selected plots into a single
  patchwork layout. If FALSE, returns a list of individual plots.

- title:

  Optional character string for the combined plot title.

- node_labels:

  Optional character vector of node names used for the heatmap
  row/column tick labels (plots 1-2). If NULL (default), labels are
  taken from the `dimnames` of the influence matrix `x$A` when present,
  otherwise `1:m`.

- period:

  Optional integer selecting a single time slice for the heatmaps when
  the fit uses dynamic influence covariates (`x$A` is `m x m x T`). If
  NULL (default), dynamic heatmaps show the time-average across all
  periods (clearly labeled as such).

- theme_base:

  A ggplot2 theme applied to all plots. Default is `theme_bw()`.

- ...:

  Additional arguments (unused).

## Value

When `combine = TRUE` and multiple plots are requested, a `patchwork`
object. When a single plot is requested, a `ggplot` object. When
`combine = FALSE`, a named list of `ggplot` objects.

## Examples

``` r
# \donttest{
dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson",
           calc_se = TRUE, seed = 1)
plot(fit)                                  # influence heatmaps + distributions

plot(fit, which = 1:6, title = "SIR Diagnostics")

plots <- plot(fit, which = c(1, 6), combine = FALSE)
plots$A_heatmap

# }
```
