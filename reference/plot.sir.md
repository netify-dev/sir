# Diagnostic Plots for a Fitted SIR Model

Produces a selection of diagnostic plots for model assessment. By
default, plots 1-4 (influence matrix heatmaps and distributions) are
shown. Use the `which` argument to select specific plots. All plots use
`ggplot2` and are combined via `patchwork` when `combine = TRUE`.

## Usage

``` r
# S3 method for class 'sir'
plot(
  x,
  which = 1:4,
  combine = TRUE,
  title = NULL,
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
if (FALSE) { # \dontrun{
model <- sir(Y, W, X, family = "poisson")

# default: influence heatmaps and distributions
plot(model)

# all plots combined with a title
plot(model, which = 1:6, title = "SIR Diagnostics")

# individual plots for custom arrangement
plots <- plot(model, which = c(1, 6), combine = FALSE)
plots$A_heatmap
plots$coef_plot
} # }
```
