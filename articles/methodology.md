# SIR Methodology

## The Social Influence Regression model

This vignette is the mathematical companion to
[`vignette("sir_overview")`](https://netify-dev.github.io/sir/articles/sir_overview.md).
Substantive readers can use the overview first, then return here for the
model scale, identification, and inference assumptions. The recurring
symbols are:

| Symbol/input             | Meaning                                          |
|:-------------------------|:-------------------------------------------------|
| `Y`                      | outcome network                                  |
| `X`                      | lagged network signal                            |
| `W`                      | covariates that build the influence matrices     |
| `Z`                      | exogenous dyadic covariates with direct effects  |
| `A`, `B`                 | sender-side and receiver-side influence matrices |
| `alpha`, `beta`, `theta` | coefficients for `W` in `A`, `W` in `B`, and `Z` |

Social Influence Regression (SIR) models directed relational data
observed over time. For a network of $`m`$ actors observed at $`T`$ time
points, let $`Y_{ijt}`$ denote the relation sent from actor $`i`$ to
actor $`j`$ at time $`t`$, with conditional mean
$`\mu_{ijt} = g^{-1}(\eta_{ijt})`$ for a link $`g`$ (log for Poisson,
identity for Normal, logit for Binomial).

The linear predictor combines exogenous direct effects with a bilinear
influence term:

``` math
\eta_{ijt} = \theta^\top z_{ijt} + (A X_t B^\top)_{ij}
    = \theta^\top z_{ijt} + \sum_{k,l} A_{ik} X_{klt} B_{jl}.
```

where

- $`z_{ijt}`$ are exogenous dyadic covariates with direct-effect
  coefficients $`\theta`$;
- $`X_t`$ is the lagged network state that carries influence (typically
  $`X_t = \log(Y_{t-1}+1)`$ for counts);
- $`A`$ (sender-side) and $`B`$ (receiver-side) are
  $`m \times m`$**influence matrices**.

The defining SIR step is to **parameterize the influence matrices by
covariates** $`W_r`$ (the slices of `W`):

``` math
A = \sum_{r=1}^{p} \alpha_r W_r,
\qquad
B = \sum_{r=1}^{p} \beta_r W_r.
```

Substituting gives the compact form
$`\eta_{ijt} = \theta^\top z_{ijt} + \alpha^\top \tilde X_{ijt} \beta`$,
where the reduced $`p \times p`$ matrix $`\tilde X_{ijt}`$ (distinct
from the $`m \times m`$ lagged state $`X_t`$) has entries
$`(\tilde X_{ijt})_{rs} = \sum_{k,l} W_{r,ik} X_{klt} W_{s,jl}`$. The
two forms are identical; the explicit double-sum above is the one used
throughout the other vignettes.

Equivalently this is a **rank-at-most-one matrix regression**: the
otherwise arbitrary $`p \times p`$ coefficient matrix $`C`$ multiplying
$`\tilde X`$ is restricted to $`C = \alpha\beta^\top`$.

In regression terms, $`\tilde X_{ijt}`$ is the table of all
sender-channel by receiver-channel lagged signals for dyad $`(i,j)`$ at
time $`t`$. SIR does not estimate an unrestricted coefficient for every
entry in that table; it restricts those interaction coefficients to
products of sender-side weights and receiver-side weights.

## Identifiability

The scale of $`\alpha`$ and $`\beta`$ is not separately identified: for
any $`c \neq 0`$, $`(c\alpha,\ \beta/c)`$ gives the identical
$`C = \alpha\beta^\top`$ and hence the same likelihood. The package
fixes $`\alpha_1 = 1`$ to resolve this, so the reported
$`\alpha_2,\dots,\alpha_p`$ are read relative to that baseline and
$`\alpha_1`$ is omitted from the coefficient table. The scale-invariant
target is the rank-at-most-one matrix $`C = \alpha\beta^\top`$ (and the
fitted means), subject to adequate design rank and signal. When
assessing recovery, compare $`C`$ rather than $`\alpha,\beta`$
separately. In applied work, put a theoretically central, nonzero
baseline channel first in `W`; all reported `alpha` coefficients are
relative to that first slice.

## Estimation

Two estimation approaches are available:

- **Alternating GLM/IRLS** (`method = "ALS"`, the default): alternates
  between updating $`(\theta, \alpha)`$ with $`\beta`$ fixed and
  $`(\theta, \beta)`$ with $`\alpha`$ fixed. Holding one influence
  vector fixed makes the model linear in the other, so each update is a
  GLM/IRLS subproblem. The `"ALS"` method name is retained for API
  compatibility.
- **optim**: maximizes the full log-likelihood jointly using BFGS with
  analytic gradients.

The alternating GLM/IRLS engine is generally more stable and is the
default.

## Standard errors

Several uncertainty summaries are available:

- **Classical**: the inverse Hessian, valid when the model is correctly
  specified, the Hessian is well conditioned, and dyad-period scores are
  independent.
- **Robust (HC0)**: a sandwich estimator for heteroskedasticity or
  overdispersion; it does not model shared-actor dyadic dependence.
- **Cluster**: a multiway cluster-robust sandwich on sender, receiver,
  and time margins; this is the usual starting point for directed
  network data when the Hessian bread is stable.
- **Bootstrap/jackknife**: `boot_sir(type = "block")` resamples whole
  periods as independent blocks, preserving within-period network
  dependence but not serial dependence, while `boot_sir(type = "dyad")`
  is a delete-one-actor jackknife.

If the Hessian is numerically singular the package falls back to a
generalized inverse and warns. A *near*-singular (but invertible)
Hessian will not trigger a warning, so also watch for a large gap
between the classical and robust SEs as a practical signal of weak
identification. Cluster SEs address dependence in the score, not weak
identification in the Hessian; when `se_reliable` is `FALSE`, prefer
refitting, simplifying the model, or using the dyad jackknife as a
sensitivity check. See
[`vignette("sir_inference")`](https://netify-dev.github.io/sir/articles/sir_inference.md)
for worked examples.

## When to use SIR

SIR is appropriate when you have a directed network observed over time
and want to explain lagged, model-implied *temporal influence* – how the
prior network predicts future ties through observed channels – using
**observed** actor/dyad covariates. It is complementary to the main
latent relational models:

- **SIR**: influence is a function of measured covariates `W` (alliance,
  distance, shared membership, …). Use it when you can name the features
  you think drive influence and want interpretable coefficients for
  them.
- **Latent space models (LSM)**: place actors in an unobserved geometric
  space; good for visualizing positions/clustering when you cannot
  specify the drivers.
- **AMEN / latent factor models**: capture additive sender/receiver
  effects and multiplicative (stochastic-equivalence) structure with
  latent factors; good for flexible dependence when covariate
  explanation is secondary.
- **ERGM**: models the network via local configuration statistics
  (reciprocity, triangles); a cross-sectional generative model, less
  geared to temporal covariate-driven influence.

A key limitation to weigh: **SIR has no contemporaneous
latent-dependence term.** Influence enters only through the *lagged*
state $`X_t`$, and the model carries no latent sender/receiver effects
or multiplicative factors for the *residual* dyadic structure within a
time period. Consequently SIR’s inference assumes dyadic independence
given the covariates — any leftover within-period dependence
(reciprocity, degree heterogeneity, stochastic equivalence) is unmodeled
and can bias the reported standard errors. AMEN and latent-space models
are preferable when that residual relational dependence is itself of
interest, or when it is strong enough that ignoring it would distort
uncertainty estimates; in applied work multiway cluster SEs and
`boot_sir(type = "dyad")` are the package’s dyad-aware sensitivity
checks. HC0 robust SEs only address heteroskedasticity/overdispersion
under dyad-period score independence.

For causal language, SIR needs the same discipline as any longitudinal
observational design: temporal ordering, exogeneity of the included
covariates, plausible control for time-varying confounding, and an
intervention or scenario that is coherent for `X`, `W`, and `Z`. Without
those assumptions, interpret the model as conditional temporal
association.

An operational noncausal report should state how `X` was constructed,
list the direct covariates in `Z`, report cluster-robust intervals (with
dyad jackknife or bootstrap checks for key coefficients), and present
named model-implied scenarios as sensitivity comparisons on the response
scale rather than as interventions.

SIR also handles **undirected** networks (set `symmetric = TRUE`, which
fixes $`B = I`$) and **bipartite / rectangular** outcomes where senders
and receivers are distinct populations. The default bipartite fit
estimates sender-side influence only; supplying `W_recv` estimates a
full bilinear sender-and-receiver model, with dyad jackknife inference
recommended; see
[`vignette("sir_extensions")`](https://netify-dev.github.io/sir/articles/sir_extensions.md).

Reach for SIR when the question is “*which observable features make
actors influential over time, and on whom?*”
