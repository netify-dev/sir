
#' Social Influence Regression (SIR) Model
#'
#' Main function to fit the SIR model. It supports Poisson, Normal, and Binary outcomes.
#' Estimation is performed using optimized C++ routines for speed.
#'
#' This implements the model described in the paper:
#' mu_{i,j,t} = theta^T z_{i,j,t} + alpha^T X_{i,j,t} beta
#'
#' @param Y (m x m x T) array of outcomes.
#' @param W (m x m x p) array of influence covariates.
#' @param X (m x m x T) array, typically lagged outcomes for the bilinear part.
#' @param Z (m x m x q x T) or (m x m x T) array of exogenous covariates.
#' @param family Distribution family ("poisson", "normal", "binomial").
#' @param method Fitting method ("ALS" or "optim"). ALS is generally recommended for stability.
#' @param calcSE Logical, whether to calculate standard errors.
#' @param ... Additional arguments passed to the fitting function (e.g., trace, tol, max_iter).
#' @return An object of class "sir" containing the estimates, standard errors, influence matrices (A, B), and fit statistics.
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


  # Prepare output
  result <- list(
    summ = summ,
    A    = A,
    B    = B,
    ll   = ll,
    NLL  = NLL,
    family = family,
    method = method,
    history = list(ALPHA=mod$ALPHA, BETA=mod$BETA, THETA=mod$THETA, DEV=mod$DEV),
    convergence = if (method=="optim") mod$convergence == 0 else (mod$iterations < 100) # Simplified convergence check
  )

  if (family == "normal") {
      result$sigma2 <- sigma2_est
  }

  class(result) <- "sir"
  return(result)
}

#' Print method for SIR objects
#' @export
print.sir <- function(x, ...) {
    cat("Social Influence Regression (SIR) Model
")
    cat("---------------------------------------
")
    cat("Family:    ", x$family, "
")
    cat("Method:    ", x$method, "
")
    cat("Log-Likelihood:", sprintf("%.4f", x$ll), "
")
    if (x$family == "normal") {
        cat("Sigma^2:   ", sprintf("%.4f", x$sigma2), "
")
    }
    cat("Converged: ", x$convergence, "
")
    cat("
Coefficients:
")
    # Print the summary table
    if ("se" %in% colnames(x$summ)) {
        # Create a printable version with rounded values
        printable_summ <- round(x$summ, 4)
        print(printable_summ)
    } else {
        print(round(x$summ, 4))
    }
}
