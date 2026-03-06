# Model framework and methodology

## The problem

Actors in political and social networks do not behave independently.
When a state initiates conflict with a particular target, that action
may reshape how other states behave toward the same or related targets
in subsequent periods. Alliance partners may follow suit, rivals may
escalate, and third parties may recalibrate their own strategies in
response. These dynamics generate higher-order dependencies across the
network that standard regression approaches cannot capture. The
fundamental challenge is not merely that observations are correlated,
but that the structure of influence itself is substantively meaningful:
understanding *who* influences *whom*, and through what channels, is
often the central question.

Traditional latent variable models for networks describe the overall
structure of interactions, positioning actors in a social space based on
transitivity or stochastic equivalence. However, while these models can
effectively characterize broad network patterns, they frequently fall
short in providing detailed explanations for the specific influence that
actors exert on one another. The factors driving influence are left
unexplored, attributed to latent dimensions rather than to the
observable actor-level and dyad-level covariates that substantive
theories emphasize. To address this limitation, the Social Influence
Regression (SIR) model regresses influence patterns directly on
observable covariates. The model operates on longitudinal network data
(a time series of $n \times n$ relational matrices) and estimates how
past interactions across the network predict current outcomes as a
function of covariates such as alliances, trade ties, and geographic
proximity.

The methodological framework is introduced in:

> Minhas, S. & Hoff, P.D. (2025). *Decomposing Network Dynamics: Social
> Influence Regression.* Political Analysis.

## Model specification

Let $Y = \{ Y_{t}:t = 1,\ldots,T\}$ be a time series of $n \times n$
relational matrices, where $y_{i,j,t}$ represents the directed outcome
from node $i$ to node $j$ at time $t$. The SIR model specifies:

$$\mu_{i,j,t} = {\mathbf{θ}}^{\top}\mathbf{z}_{i,j,t} + {\mathbf{α}}^{\top}{\widetilde{X}}_{i,j,t}{\mathbf{β}}$$

The first term is a standard regression: $\mathbf{z}_{i,j,t}$ collects
exogenous covariates for dyad $(i,j)$ at time $t$ (geographic distance,
alliance status, trade flows), and $\mathbf{θ}$ gives their direct
effects on the outcome. The second term is where the model departs from
standard approaches. The quantity ${\widetilde{X}}_{i,j,t}$ is
constructed from the lagged network state interacted with influence
covariates $\mathbf{W}$, and the bilinear form
${\mathbf{α}}^{\top}{\widetilde{X}}_{i,j,t}{\mathbf{β}}$ allows each
dyad’s outcome to depend on the *entire* prior network, weighted by
sender ($\mathbf{α}$) and receiver ($\mathbf{β}$) influence parameters.
This structure is what allows the model to capture third-order
dependencies: how one actor’s past behavior toward a third party
predicts another actor’s current behavior.

For any pair of actors $(i,i\prime)$, the entry $a_{i,i\prime}$ tells us
how predictive $i\prime$’s past sending behavior is of $i$’s current
sending behavior. If $a_{\text{GBR},\text{USA}} > 0$ in a conflict
network, it indicates that countries the USA initiated conflict with in
period $t - 1$ tend to also face conflict from the UK in period $t$. The
influence is directional and asymmetric: the USA’s conflict behavior is
predictive of the UK’s, but the reverse need not hold with the same
magnitude.

### Influence matrices

The influence parameters $a_{i,i\prime}$ and $b_{j,j\prime}$ are not
estimated freely for every pair, as that would require
$O\left( n^{2} \right)$ parameters. Instead, the model explains
influence in terms of covariates:

$$a_{i,i\prime} = \alpha^{\top}w_{i,i\prime}\qquad b_{j,j\prime} = \beta^{\top}w_{j,j\prime}$$

where $w_{i,i\prime}$ is a vector of covariates describing the
relationship between actors $i$ and $i\prime$ (distance, alliance
status, trade ties). The full influence matrices are then:

$$\mathbf{A} = \sum\limits_{r = 1}^{p}\alpha_{r}\mathbf{W}_{r}\qquad\mathbf{B} = \sum\limits_{r = 1}^{p}\beta_{r}\mathbf{W}_{r}$$

with $\alpha_{1} = 1$ fixed for identifiability. This parameterization
brings the parameter count down to $O(p)$ where $p$ is the number of
influence covariates, typically a handful. The $\alpha$ and $\beta$
coefficients tell us which covariates matter for influence and by how
much: a positive $\alpha$ on the alliance covariate, for instance, would
indicate that allied countries tend to initiate conflict with the same
targets. This emphasis on covariate-driven explanation is what
distinguishes the SIR framework from latent variable approaches. Rather
than describing influence through unobserved dimensions, the model links
influence directly to measured actor and dyad attributes.

### Distribution families

The framework is based on a generalized bilinear model and extends
naturally to different outcome types:

| Family   | Link                         | Example                                         |
|:---------|:-----------------------------|:------------------------------------------------|
| Poisson  | $g(\mu) = \log(\mu)$         | Monthly conflict event counts between countries |
| Normal   | $g(\mu) = \mu$               | Bilateral trade volumes                         |
| Binomial | $g(\mu) = \text{logit}(\mu)$ | Presence or absence of a diplomatic tie         |

### Identifiability

The bilinear term ${\mathbf{α}}^{\top}X{\mathbf{β}}$ has a scale
ambiguity: multiplying $\mathbf{α}$ by $1/c$ and $\mathbf{β}$ by $c$
yields the same product for any nonzero scalar $c$. To resolve this, the
first element of $\mathbf{α}$ is fixed at 1. The package handles this
constraint automatically during estimation.

A simpler alternative is to fix $\mathbf{B} = \mathbf{I}$ entirely
(`fix_receiver = TRUE`), which removes the identification issue and
reduces the model to a standard GLM. This is a useful starting point
when the research question concerns sender-side influence alone, and it
produces well-conditioned standard errors without requiring bootstrap
corrections.

## Estimation

Estimating $\{{\mathbf{θ}},{\mathbf{α}},{\mathbf{β}}\}$ jointly is
difficult because the model is bilinear in $\mathbf{α}$ and
$\mathbf{β}$. The package addresses this with an iterative block
coordinate descent algorithm that exploits a key structural property:
for fixed $\mathbf{β}$, the model is linear in $\mathbf{θ}$ and
$\mathbf{α}$ (and vice versa). The procedure initializes $\mathbf{β}$,
then alternates between two steps: first, fixing $\mathbf{β}$ and
estimating $({\mathbf{θ}},{\mathbf{α}})$ via GLM; second, fixing
$\mathbf{α}$ and estimating $({\mathbf{θ}},{\mathbf{β}})$ via GLM.
Iteration continues until the relative change in deviance falls below a
tolerance. Each sub-problem is a standard generalized linear model
solved by iterative weighted least squares, so the full estimation
reduces to a sequence of low-dimensional optimizations. This is
substantially faster than the Bayesian approach originally used for
bilinear network autoregressions.

The package also provides a direct BFGS method (`method = "optim"`) that
optimizes all parameters simultaneously using analytical gradients
computed via C++. This approach can converge faster for small networks
but tends to be less stable when $p$ is large.

## Inference

Standard errors come from the Hessian of the log-likelihood at the MLE.
The package computes classical standard errors from the observed
information matrix $H^{- 1}$ as well as robust (sandwich) standard
errors $H^{- 1}SH^{- 1}$, where $S$ is the empirical score covariance.
The sandwich estimator remains valid under model misspecification.

The Hessian can be ill-conditioned in bilinear models, and the package
warns when this occurs. In such cases, bootstrap standard errors via
[`boot_sir()`](https://netify-dev.github.io/sir/reference/boot_sir.md)
provide a more reliable basis for inference. The bootstrap supports both
block resampling of time periods and parametric simulation from the
fitted model.

## Choosing a model configuration

| Consideration                                      | Full bilinear | Fixed receiver |
|:---------------------------------------------------|:-------------:|:--------------:|
| Both sender and receiver influence                 |       ✓       |                |
| Standard GLM with well-conditioned standard errors |               |       ✓        |
| Bilinear identification resolved                   |               |       ✓        |
| Richer influence structure                         |       ✓       |                |
| Fewer parameters                                   |               |       ✓        |

For exploratory analysis, the fixed-receiver specification is a
reasonable starting point: fewer parameters, no identification issues,
and proper standard errors from the GLM. The full bilinear model is
appropriate when there is reason to believe that both sender and
receiver channels contribute to influence dynamics.

| Consideration                       | ALS | Direct optimization |
|:------------------------------------|:---:|:-------------------:|
| High-dimensional problems           |  ✓  |                     |
| Small problems (fast convergence)   |     |          ✓          |
| Numerical stability                 |  ✓  |                     |
| Simultaneous parameter optimization |     |          ✓          |

## Citation

If you use this package in your research, please cite:

> Minhas, S. & Hoff, P.D. (2025). Decomposing Network Dynamics: Social
> Influence Regression. Political Analysis.

The package source code and documentation are available at
<https://github.com/netify-dev/sir>.
