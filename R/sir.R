
#' Social Influence Regression (SIR) Model
#'
#' @description
#' Fits a Social Influence Regression model for network data with social influence effects.
#' The SIR model captures how network connections influence outcomes through bilinear 
#' interaction terms, allowing for both sender and receiver effects in directed networks.
#' 
#' The model decomposes network influence into two components:
#' \itemize{
#'   \item \strong{Sender effects (A matrix)}: How much influence node i exerts on others
#'   \item \strong{Receiver effects (B matrix)}: How susceptible node j is to influence from others
#' }
#'
#' @details
#' The SIR model specifies the expected outcome for the directed edge from node i to node j 
#' at time t as:
#' 
#' \deqn{\mu_{i,j,t} = \theta^T z_{i,j,t} + \sum_{k,l} X_{k,l,t} A_{i,k} B_{j,l}}
#' 
#' Where:
#' \itemize{
#'   \item \eqn{\mu_{i,j,t}} is the expected value of the outcome Y_ijt
#'   \item \eqn{\theta} is a q-dimensional vector of coefficients for exogenous covariates
#'   \item \eqn{z_{i,j,t}} is a q-dimensional vector of exogenous covariates
#'   \item \eqn{X_{k,l,t}} represents the network state (often lagged Y) that carries influence
#'   \item \eqn{A_{i,k}} represents how node i is influenced by the behavior of node k
#'   \item \eqn{B_{j,l}} represents how node j's reception is affected by node l's position
#' }
#' 
#' The bilinear term \eqn{\sum_{k,l} X_{k,l,t} A_{i,k} B_{j,l}} captures network influence
#' and can be parameterized using influence covariates W through:
#' \itemize{
#'   \item \eqn{A = \alpha_0 I + \sum_{r=1}^{p} \alpha_r W_r} (sender effects)
#'   \item \eqn{B = \beta_0 I + \sum_{r=1}^{p} \beta_r W_r} (receiver effects)
#' }
#' 
#' This parameterization reduces the number of parameters from O(m²) to O(p), where p << m.
#' 
#' @section Estimation Methods:
#' 
#' \strong{Alternating Least Squares (ALS):}
#' \itemize{
#'   \item Iteratively optimizes A given B, then B given A
#'   \item Generally more stable for high-dimensional problems
#'   \item Better for sparse networks or when p is large
#'   \item May converge to local optima
#' }
#' 
#' \strong{Direct Optimization (optim):}
#' \itemize{
#'   \item Uses BFGS to optimize all parameters simultaneously
#'   \item Can be faster for small problems
#'   \item May provide better solutions when good starting values are available
#'   \item More prone to numerical issues in high dimensions
#' }
#'
#' @section Distribution Families:
#' 
#' \strong{Poisson:} For count data (e.g., number of interactions)
#' \itemize{
#'   \item Link function: log
#'   \item Variance function: V(μ) = μ
#'   \item Use when: Y_ijt represents counts
#' }
#' 
#' \strong{Normal:} For continuous data (e.g., trade volumes, distances)
#' \itemize{
#'   \item Link function: identity
#'   \item Variance function: V(μ) = σ²
#'   \item Use when: Y_ijt is continuous and approximately normal
#' }
#' 
#' \strong{Binomial:} For binary data (e.g., presence/absence of ties)
#' \itemize{
#'   \item Link function: logit
#'   \item Variance function: V(μ) = μ(1-μ)
#'   \item Use when: Y_ijt is binary (0/1)
#' }
#'
#' @param Y A three-dimensional array of dimensions (m x m x T) containing the network 
#'   outcomes. Y[i,j,t] represents the directed outcome from node i to node j at time t.
#'   Can contain NA values for missing observations. The diagonal (self-loops) can be 
#'   included or excluded depending on the application.
#'   
#' @param W Optional three-dimensional array of dimensions (m x m x p) containing 
#'   influence covariates used to parameterize the A and B matrices. W[i,j,r] represents
#'   the r-th influence covariate for the edge from i to j. Common choices include:
#'   \itemize{
#'     \item Graph Laplacians or adjacency matrices from other networks
#'     \item Geographic or social distance matrices
#'     \item Node-level covariates expanded to edge-level
#'   }
#'   If NULL or p=0, the model uses only identity matrices (no network influence structure).
#'   
#' @param X Optional three-dimensional array of dimensions (m x m x T) representing the 
#'   network state that carries influence. Typically this is a lagged version of Y 
#'   (e.g., X[,,t] = Y[,,t-1]). If NULL and W is provided, an error is thrown.
#'   X determines which network patterns influence future outcomes.
#'   
#' @param Z Optional array of exogenous covariates. Can be either:
#'   \itemize{
#'     \item 3D array (m x m x T): Single covariate varying across edges and time
#'     \item 4D array (m x m x q x T): Multiple (q) covariates
#'   }
#'   Examples include dyadic covariates (trade agreements, geographic distance) or
#'   node-level attributes (GDP, population) expanded to edge-level.
#'   
#' @param family Character string specifying the distribution family and link function.
#'   Must be one of "poisson", "normal", or "binomial". The choice depends on the
#'   nature of your outcome variable.
#'   
#' @param method Character string specifying the estimation method. Either "ALS" 
#'   (Alternating Least Squares) or "optim" (direct optimization via BFGS).
#'   Default is "ALS" which is generally more stable.
#'   
#' @param calcSE Logical indicating whether to calculate standard errors for the
#'   parameters. Standard errors are computed using the observed information matrix.
#'   Setting to FALSE speeds up computation when uncertainty quantification is not needed.
#'   
#' @param ... Additional arguments passed to the fitting functions:
#'   \itemize{
#'     \item \code{trace}: Logical or integer controlling output verbosity
#'     \item \code{tol}: Convergence tolerance (default 1e-6)
#'     \item \code{max_iter}: Maximum iterations (default 100 for ALS, 1000 for optim)
#'     \item \code{init_method}: Initialization method ("smart" or "random")
#'   }
#'   
#' @return An object of class "sir" containing:
#'   \item{A}{The m x m sender effects matrix}
#'   \item{B}{The m x m receiver effects matrix}
#'   \item{theta}{Coefficients for exogenous covariates Z}
#'   \item{alpha}{Coefficients for sender influence covariates W}
#'   \item{beta}{Coefficients for receiver influence covariates W}
#'   \item{summ}{Summary table with coefficients and standard errors}
#'   \item{ll}{Log-likelihood at convergence}
#'   \item{AIC}{Akaike Information Criterion}
#'   \item{BIC}{Bayesian Information Criterion}
#'   \item{convergence}{Logical indicating successful convergence}
#'   \item{iterations}{Number of iterations until convergence}
#'   \item{history}{Iteration history of parameters and deviance (if trace=TRUE)}
#'   \item{call}{The matched call}
#'   
#' @references
#' Minhas, S., Hoff, P. D., & Ward, M. D. (2019). Inferential approaches for 
#' network analysis: AMEN for latent factor models. Political Analysis, 27(2), 208-222.
#' 
#' Hoff, P. D. (2021). Additive and multiplicative effects network models. 
#' Statistical Science, 36(1), 34-50.
#' 
#' @examples
#' \dontrun{
#' # Generate example network data
#' set.seed(123)
#' m <- 20  # number of nodes
#' T <- 10  # number of time periods
#' 
#' # Create synthetic network with influence
#' Y <- array(rpois(m*m*T, lambda=2), dim=c(m,m,T))
#' 
#' # Use lagged Y as influence carrier
#' X <- array(0, dim=c(m,m,T))
#' X[,,2:T] <- Y[,,1:(T-1)]
#' 
#' # Create influence covariates (e.g., based on node attributes)
#' W <- array(rnorm(m*m*2), dim=c(m,m,2))
#' 
#' # Fit SIR model with Poisson family
#' model <- sir(Y=Y, W=W, X=X, family="poisson", method="ALS")
#' 
#' # View results
#' summary(model)
#' plot(model)
#' 
#' # Extract influence matrices
#' A_matrix <- model$A  # Sender effects
#' B_matrix <- model$B  # Receiver effects
#' }
#' 
#' @importFrom MASS ginv
#' @export
sir <- function(Y, W=NULL, X=NULL, Z=NULL, family, method="ALS", calcSE=TRUE, ...) {

  if (!family %in% c("poisson", "normal", "binomial")) {
      stop("Invalid family specified.")
  }

  # Input validation and dimension checks
  dimsY <- dim(Y)
  m <- dimsY[1]
  T_len <- dimsY[3]

  if (length(dimsY) != 3 || dimsY[1] != dimsY[2]) {
      stop("Y must be a 3D array (m x m x T).")
  }

  if (!is.null(W)) {
      p <- dim(W)[3]
      if (is.null(X)) stop("X must be provided if W is provided.")
      if (any(dimsY[1:2] != dim(W)[1:2]) || any(dimsY != dim(X))) {
          stop("Dimensions of Y (m x m x T), W (m x m x p), and X (m x m x T) must align.")
      }
  } else {
      p <- 0
      # Create dummy empty arrays for Cpp functions if p=0.
      W <- array(0, dim=c(m, m, 0))
      # X is still needed if p=0 for eta_tab calculation, ensure it matches Y dims.
      if (is.null(X)) X <- array(0, dim=c(m, m, T_len))
  }

  # Handle Z dimensions (Allow 3D or 4D)
  if (!is.null(Z)) {
      dimsZ <- dim(Z)
      if (length(dimsZ) == 3) {
          if (all(dimsZ == dimsY)) {
              # Assuming (m x m x T) and q=1
              Z <- array(Z, dim=c(m, m, 1, T_len))
              q <- 1
          } else {
              stop("3D Z array must have dimensions (m x m x T) matching Y.")
          }
      } else if (length(dimsZ) == 4) {
          if (all(dimsZ[c(1,2,4)] == dimsY)) {
              q <- dimsZ[3]
          } else {
              stop("4D Z array must have dimensions (m x m x q x T) matching Y.")
          }
      } else {
          stop("Z must be a 3D or 4D array.")
      }
  } else {
      q <- 0
      # Z remains NULL for function calls that check for it
  }


  # ---- Fitting ----
  if (method == "ALS") {
    mod <- sir_alsfit(Y, W, X, Z, family, ...)
  } else if (method == "optim") {
    mod <- sir_optfit(Y, W, X, Z, family, ...)
  } else {
    stop("Invalid method specified. Choose 'ALS' or 'optim'.")
  }

  tab <- mod$tab

  # ---- Post-processing ----

  # Reconstruct final alpha, beta
  if (q > 0) {
      theta <- tab[1:q]
  } else {
      theta <- numeric(0)
  }

  if (p > 0) {
      if (p > 1) {
          a_start = q + 1
          a_end = q + p - 1
          a <- tab[a_start:a_end]
          alpha <- c(1, a)
      } else {
          a <- numeric(0)
          alpha <- c(1)
      }
      b_start = q + p
      beta <- tab[b_start:length(tab)]
  } else {
      alpha <- numeric(0)
      beta <- numeric(0)
  }


  # Influence matrices A, B (using Cpp optimized functions)
  if (p > 0) {
      A <- cpp_amprod_W_v(W, alpha)
      B <- cpp_amprod_W_v(W, beta)

      # Optional normalization (as in original scripts)
      # Check mean(diag(A)) safely
      diag_A_mean <- mean(diag(A), na.rm=TRUE)
      if (!is.na(diag_A_mean) && diag_A_mean != 0) {
        A <- A * sign(diag_A_mean)
      }
      diag_B_mean <- mean(diag(B), na.rm=TRUE)
      if (!is.na(diag_B_mean) && diag_B_mean != 0) {
        B <- B * sign(diag_B_mean)
      }
      diag(A) <- 0
      diag(B) <- 0

      if (!is.null(rownames(Y))) {
        rownames(A) <- colnames(A) <- rownames(Y)
        rownames(B) <- colnames(B) <- rownames(Y)
      }
  } else {
      A <- matrix(0, m, m)
      B <- matrix(0, m, m)
  }


  # Log-Likelihood and Variance Estimation
  # mll_sir calculates NLL (assuming sigma=1 for normal).
  NLL <- mll_sir(tab, Y, W, X, Z, family)
  ll <- -NLL
  sigma2_est <- 1.0 # Default

  if (family == "normal") {
      # Recalculate LL and estimate sigma^2 for Gaussian model
      N_obs <- sum(!is.na(Y))
      df <- N_obs - length(tab)
      # RSS is stored in DEV[last, 2] for ALS, or 2*NLL for optim (if sigma=1)
      if (method == "ALS") {
          RSS <- mod$DEV[nrow(mod$DEV), 2]
      } else {
          # If optimized assuming sigma=1, NLL = 0.5 * RSS.
          RSS <- 2 * NLL
      }
      sigma2_est <- RSS / df
      # Full Gaussian LL
      ll <- -0.5 * N_obs * log(2 * pi * sigma2_est) - 0.5 * RSS / sigma2_est
  }


  # ---- Standard Errors ----
  summ <- data.frame(coef=tab)
  if(calcSE){
    # Use the Rcpp backend for Hessian calculation (Hessian of NLL)
    Z_list <- prepare_Z_list(Z)
    gH <- cpp_mll_gH(tab, Y, W, X, Z_list, family)

    H <- gH$hess

    # Adjust Hessian/Score for Gaussian sigma^2
    if (family == "normal") {
        # H_true = H_cpp / sigma^2
        H <- H / sigma2_est
        # S_true = S_cpp / sigma^4 (if S_cpp calculated using (mu-y) residuals)
        # Note: cpp_mll_gH calculates S using d(NLL)/d(param), which already includes the sigma scaling if NLL definition included it.
        # However, the Cpp implementation assumes sigma=1. We adjust here.
        S <- gH$shess / (sigma2_est^2)
    } else {
        S <- gH$shess
    }


    # Check if Hessian is invertible (should be positive definite for NLL minimum)
    H_inv <- tryCatch(solve(H), error = function(e) {
        warning("Hessian matrix is singular at the optimum, trying generalized inverse (MASS::ginv). Standard errors might be unreliable.")
        tryCatch(MASS::ginv(H), error = function(e2) {
            warning("Generalized inverse failed. Cannot compute standard errors.")
            return(NULL)
        })
    })

    if (!is.null(H_inv)) {
        # Classical SE
        se_diag <- diag(H_inv)
        # Ensure non-negative variance estimates
        se_diag[se_diag < 0] <- NA
        se <- sqrt(se_diag)

        # Robust (Sandwich) SE
        # Sandwich = H_inv %*% S %*% H_inv
        rse <- tryCatch({
            rse_diag <- diag(H_inv %*% S %*% H_inv)
            rse_diag[rse_diag < 0] <- NA
            sqrt(rse_diag)
        }, error = function(e) {
            warning("Sandwich estimator calculation failed.")
            return(rep(NA, length(se)))
        })

        summ$se <- se
        summ$rse <- rse
        summ$t_se <- summ$coef / summ$se
        summ$t_rse <- summ$coef / summ$rse
    } else {
        summ$se <- NA
        summ$rse <- NA
        summ$t_se <- NA
        summ$t_rse <- NA
    }
  }

  # Assign meaningful row names
  Z_names <- if (q > 0 && !is.null(dimnames(Z)[[3]])) dimnames(Z)[[3]] else paste0("Z", 1:q)
  W_names <- if (p > 0 && !is.null(dimnames(W)[[3]])) dimnames(W)[[3]] else paste0("W", 1:p)

  theta_names <- if (q>0) paste0("(Z) ", Z_names) else NULL
  alpha_names <- if (p>1) paste0("(alphaW) ", W_names[-1]) else NULL
  beta_names  <- if (p>0) paste0("(betaW) ", W_names) else NULL

  rownames(summ) <- c(theta_names, alpha_names, beta_names)


  # Count observations
  N_obs <- sum(!is.na(Y))
  
  # Prepare output
  result <- list(
    summ = summ,
    A    = A,
    B    = B,
    ll   = ll,
    NLL  = NLL,
    family = family,
    method = method,
    tab = tab,  # Store full parameter vector
    theta = theta,
    alpha = alpha,
    beta = beta,
    p = p,  # Number of influence covariates
    q = q,  # Number of exogenous covariates
    T = T_len,  # Number of time periods
    nobs = N_obs,  # Number of observations
    iterations = if (!is.null(mod$iterations)) mod$iterations else NA,
    history = list(ALPHA=mod$ALPHA, BETA=mod$BETA, THETA=mod$THETA, DEV=mod$DEV),
    convergence = if (method=="optim") mod$convergence == 0 else (mod$iterations < 100), # Simplified convergence check
    call = match.call()  # Store the call
  )

  if (family == "normal") {
      result$sigma2 <- sigma2_est
  }

  class(result) <- "sir"
  return(result)
}

