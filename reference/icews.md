# ICEWS Inter-State Material Conflict (Monthly)

A longitudinal directed network of monthly material-conflict event
counts between countries, derived from the Integrated Crisis Early
Warning System (ICEWS) and used to illustrate the Social Influence
Regression model of Minhas & Hoff (2025). The data are subset to the 50
countries most involved in material conflict over the sample window
(February 2005 through December 2012) to keep the bundled object small.
See the installed provenance note
`system.file("icews-provenance.md", package = "sir")` for derivation
details and repository-file checksums.

## Usage

``` r
icews
```

## Format

A named `list` with the array inputs
[`sir`](https://netify-dev.github.io/sir/reference/sir.md) expects:

- Y:

  Integer array `50 x 50 x 95`. `Y[i, j, t]` is the number of
  material-conflict events initiated by country `i` toward country `j`
  in month `t`. Rows are senders, columns receivers; the diagonal is
  stored as zero in the bundled object and ignored by
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md) for
  one-mode fits.

- X:

  Numeric array `50 x 50 x 95`, the influence-carrying network state
  `log(Y_{ij,t-1} + 1)` (the logged, lagged outcome). This is what the
  estimated influence matrices act on. The first bundled month uses a
  retained January 2005 lag from the source archive, so it cannot be
  reconstructed from the bundled `Y` alone.

- W:

  Numeric array `50 x 50 x 4` of static (time-averaged) influence
  covariates that parameterize the sender/receiver influence matrices A
  and B. Slices: `int` (intercept anchor for the \\\alpha_1 = 1\\
  identifiability constraint), `ally` (alliance ties), `verbCoop`
  (logged verbal cooperation), and `minDistLog` (log minimum geographic
  distance).

- Z:

  Numeric array `50 x 50 x 5 x 95` of exogenous dyadic covariates with
  direct effects (theta). Slices: `mConf` (lagged material conflict \\i
  \to j\\), `mConf_ji` (lagged reciprocal conflict \\j \to i\\),
  `minDistLog`, `ally`, and `verbCoop` (logged).

- countries:

  Character vector of the 50 country names (the row/column labels of
  `Y`).

- dates:

  Character vector of 95 month labels from `"2005-02-01"` through
  `"2012-12-01"`.

- metadata:

  List with source DOI, replication DOI, raw archive path, access date,
  country-selection rule, date window, and transformation notes.

## Source

Integrated Crisis Early Warning System (ICEWS), via
`replArchive/data/socRegData.rda` from the replication archive for
Minhas, S. & Hoff, P. D. (2025), "Decomposing Network Influence: Social
Influence Regression", *Political Analysis*. Replication archive DOI:
[doi:10.7910/DVN/VTFDX6](https://doi.org/10.7910/DVN/VTFDX6) . See also
the bundled provenance note:
`system.file("icews-provenance.md", package = "sir")`.

## Details

This bundled 50-country by 95-month slice is provided purely to
demonstrate the package on realistic data. It is **illustrative only**
and is *not* the estimation sample or the published results of Minhas &
Hoff (2025); fits on this slice will not reproduce the numbers reported
in the paper, and should not be cited as such.

This is the canonical applied example for the package. A typical fit:


      data(icews)
      fit <- sir(icews$Y, W = icews$W, X = icews$X, Z = icews$Z,
                 family = "poisson", seed = 1)
      summary(fit)

The covariate construction (logging verbal cooperation, the `int`
intercept slice as the influence anchor) mirrors the replication archive
for the published analysis.

## References

Minhas, S. & Hoff, P. D. (2025). Decomposing Network Influence: Social
Influence Regression. *Political Analysis*.
[doi:10.1017/pan.2025.10013](https://doi.org/10.1017/pan.2025.10013) .

## Examples

``` r
data(icews)
dim(icews$Y)            # 50 x 50 x 95
#> [1] 50 50 95
head(icews$countries)
#> [1] "UNITED STATES" "AFGHANISTAN"   "ISRAEL"        "PAKISTAN"     
#> [5] "LEBANON"       "IRAQ"         
if (FALSE) { # \dontrun{
fit <- sir(icews$Y, W = icews$W, X = icews$X, Z = icews$Z,
           family = "poisson", seed = 1)
summary(fit)
} # }
```
