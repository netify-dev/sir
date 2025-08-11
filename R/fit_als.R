
#' Fit SIR model via Alternating Least Squares (Block Coordinate Descent)
#'
#' Iteratively updates (theta, alpha) and (theta, beta) using GLMs.
#' This implementation uses optimized C++ functions to construct the design matrices.
#'
#' @param Y (m x m x T) outcome array.
#' @param W (m x m x p) influence covariates.
#' @param X (m x m x T) bilinear covariates (e.g., lagged Y).
#' @param Z (m x m x q x T) exogenous covariates.
#' @param family Distribution family ("poisson", "normal", "binomial").
#' @param trace Logical, whether to print iteration details.
#' @param tol Convergence tolerance.
#' @param max_iter Maximum number of iterations.
#' @return A list containing the estimated parameters and iteration history.
#' @importFrom stats glm lm poisson binomial coef deviance formula
#' @importFrom speedglm speedglm
#' @export
sir_alsfit <- function(Y, W, X, Z, family, trace=FALSE, tol=1e-8, max_iter=100) {
  p <- if (is.null(W)) 0 else dim(W)[3]
  q <- if (is.null(Z)) 0 else dim(Z)[3]
  m <- dim(Y)[1]
  T_len <- dim(Y)[3]

  # Determine GLM function and family object
  if (family == "normal") {
    glm_fun <- stats::lm
    family_obj <- NULL
    use_speedglm <- FALSE
  } else {
    # Use speedglm for faster fitting if available
    if (requireNamespace("speedglm", quietly = TRUE)) {
        glm_fun <- speedglm::speedglm
        use_speedglm <- TRUE
    } else {
        glm_fun <- stats::glm
        use_speedglm <- FALSE
    }
    family_obj <- switch(family,
                         "poisson" = stats::poisson(),
                         "binomial" = stats::binomial(),
                         stop("Unsupported family for GLM."))
  }

  # ---- Initialization ----
  Y_flat <- flatten_Y(Y)
  Z_flat <- flatten_Z(Z)

  if (q > 0) {
    glmData <- data.frame(Y = Y_flat, Z_flat)
    form <- formula(paste0("Y ~ -1 + ", paste(colnames(Z_flat), collapse=" + ")))

    if (family == "normal") {
        fit0 <- glm_fun(form, data=glmData)
    } else {
        fit0 <- glm_fun(form, data=glmData, family=family_obj)
    }
    theta <- coef(fit0)
    # Handle potential NA coefficients from initialization
    theta[is.na(theta)] <- 0
    dev_new <- if (family == "normal") sum(fit0$residuals^2) else deviance(fit0)

  } else {
      theta <- numeric(0)
      # Calculate null deviance if q=0 (simplified initialization)
      mean_Y <- mean(Y_flat, na.rm=TRUE)
      if (family == "normal") {
          dev_new <- sum((Y_flat - mean_Y)^2, na.rm=TRUE)
      } else if (family == "poisson") {
          # Avoid log(0) issues
          if (mean_Y <= 0) mean_Y <- 1e-6
          dev_new = 2 * sum(ifelse(Y_flat > 0, Y_flat * log(Y_flat/mean_Y), 0) - (Y_flat - mean_Y), na.rm=TRUE)
      } else if (family == "binomial") {
          p0 = mean_Y
          if (p0 <= 0) p0 = 1e-6
          if (p0 >= 1) p0 = 1 - 1e-6
          dev_new = -2 * sum(Y_flat * log(p0) + (1-Y_flat) * log(1-p0), na.rm=TRUE)
      }
  }


  # Initialize alpha, beta (small random values)
  if (p > 0) {
      # Initialize slightly larger than 1e-3 for stability
      alpha <- rnorm(p, sd=0.05)
      beta  <- rnorm(p, sd=0.05)
  } else {
      alpha <- numeric(0)
      beta <- numeric(0)
  }


  # Track iteration
  dev_old <- Inf
  THETA <- matrix(theta, nrow=1)
  ALPHA <- matrix(alpha, nrow=1)
  BETA  <- matrix(beta, nrow=1)
  DEV   <- matrix(c(dev_old, dev_new), nrow=1)

  # ---- Iterative Updates ----
  for (iter in 1:max_iter) {

    if (p == 0) break; # Stop if no bilinear part

    # 1) Fix beta, update (theta, alpha)

    # Construct design matrix Wbeta using optimized C++ function
    # This returns the flattened (m*m*T) x p matrix directly.
    Wbeta_flat <- cpp_construct_Wbeta_design(W, X, beta)
    colnames(Wbeta_flat) <- paste0("WB", 1:p)

    # Combine data and define formula
    if (q > 0) {
        X_design <- cbind(Z_flat, Wbeta_flat)
    } else {
        X_design <- Wbeta_flat
    }
    glmData <- data.frame(Y = Y_flat, X_design)
    form1 <- formula(paste0("Y ~ -1 + ", paste(colnames(X_design), collapse=" + ")))

    # Fit GLM
    start_vals <- c(theta, alpha)
    fitA <- fit_glm_wrapper(form1, glmData, family, family_obj, glm_fun, use_speedglm, start_vals, trace)

    coA  <- coef(fitA)
    # Handle potential NA coefficients during iteration
    coA[is.na(coA)] <- 0

    if (q > 0) {
        theta <- coA[1:q]
        alpha <- coA[-(1:q)]
    } else {
        alpha <- coA
    }

    # 2) Fix alpha, update (theta, beta)

    # Construct design matrix Walpha using optimized C++ function
    Walpha_flat <- cpp_construct_Walpha_design(W, X, alpha)
    colnames(Walpha_flat) <- paste0("WA", 1:p)

    # Combine data and define formula
    if (q > 0) {
        X_design <- cbind(Z_flat, Walpha_flat)
    } else {
        X_design <- Walpha_flat
    }
    glmData <- data.frame(Y = Y_flat, X_design)
    form2 <- formula(paste0("Y ~ -1 + ", paste(colnames(X_design), collapse=" + ")))

    # Fit GLM
    start_vals <- c(theta, beta)
    fitB <- fit_glm_wrapper(form2, glmData, family, family_obj, glm_fun, use_speedglm, start_vals, trace)

    coB  <- coef(fitB)
    coB[is.na(coB)] <- 0

    if (q > 0) {
        theta <- coB[1:q]
        beta  <- coB[-(1:q)]
    } else {
        beta <- coB
    }

    # Check deviance
    dev_old <- dev_new
    dev_new <- if (family == "normal") sum(fitB$residuals^2) else deviance(fitB)

    # Handle potential NA deviance
    if (is.na(dev_new)) {
        if (trace) warning("Deviance became NA, stopping iteration.")
        break
    }

    # Track history
    THETA <- rbind(THETA, theta)
    ALPHA <- rbind(ALPHA, alpha)
    BETA  <- rbind(BETA, beta)
    DEV   <- rbind(DEV, c(dev_old, dev_new))

    if(trace){
      cat("Iteration:", iter, "Dev:", sprintf("%.4f", dev_new), "Change:", sprintf("%.6f", abs(dev_old - dev_new)/(abs(dev_old) + 0.1)), "
")
    }

    # Stopping criterion
    if( abs(dev_old - dev_new) / (abs(dev_old) + 0.1) < tol ) {
        if (trace) cat("Converged.
")
        break
    }
  }

  if (iter == max_iter && p > 0) {
      warning("ALS did not converge within the maximum number of iterations.")
  }

  # Impose alpha_1=1 constraint and rescale
  if (p > 0) {
      if (abs(alpha[1]) < 1e-8) {
          warning("Estimated alpha[1] is near zero. Rescaling might be unstable. Using pseudoinverse stabilization.")
          # Stabilization if alpha[1] is too small
          a1 <- sign(alpha[1])
          if (a1 == 0) a1 <- 1 # If alpha[1] was exactly 0, pick a sign
          a1 <- a1 * max(abs(alpha[1]), 1e-6)
      } else {
          a1 <- alpha[1]
      }

      if (p > 1) {
        a <- alpha[-1] / a1
      } else {
        a <- numeric(0)
      }
      b <- beta * a1
  } else {
      a <- numeric(0)
      b <- numeric(0)
  }


  list(
    theta=theta,
    a=a,
    b=b,
    tab=c(theta, a, b),
    ALPHA=ALPHA,
    BETA=BETA,
    THETA=THETA,
    DEV=DEV,
    iterations = iter
  )
}

# Helper function to handle GLM fitting consistently and robustly
fit_glm_wrapper <- function(form, data, family, family_obj, glm_fun, use_speedglm, start_vals, trace) {
    if (family == "normal") {
        fit <- glm_fun(form, data=data)
    } else {
        # Ensure start values match the number of predictors
        if (length(start_vals) != ncol(data)-1) {
             if (trace) warning("Start values length mismatch. Using default initialization.")
             fit <- glm_fun(form, data=data, family=family_obj)
        } else {
            # Try using start values, robustly handle potential failures (common in IWLS)
            tryCatch({
                if (use_speedglm) {
                    # speedglm specific arguments
                    fit <- glm_fun(form, data=data, family=family_obj, start=start_vals)
                } else {
                    # stats::glm arguments
                    fit <- glm_fun(form, data=data, family=family_obj, start=start_vals)
                }
            }, error = function(e) {
                if (trace) warning("GLM fitting failed with start values (", e$message, "), falling back to default initialization.")
                # Fallback if starting values cause issues
                tryCatch({
                     fit <<- glm_fun(form, data=data, family=family_obj)
                }, error = function(e2) {
                    # Final fallback using stats::glm if speedglm fails entirely
                    if (use_speedglm) {
                        if (trace) warning("speedglm failed entirely, falling back to stats::glm.")
                        fit <<- stats::glm(form, data=data, family=family_obj)
                    } else {
                        stop("stats::glm failed during ALS iteration: ", e2$message)
                    }
                })
            })
        }
    }
    return(fit)
}
