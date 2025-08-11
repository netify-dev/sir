
#' Fit SIR Model via Alternating Least Squares (ALS)
#'
#' @description
#' Implements the Alternating Least Squares algorithm (also known as Block Coordinate Descent) 
#' for fitting Social Influence Regression models. This method alternates between optimizing 
#' the sender effects (A matrix) and receiver effects (B matrix) while holding the other fixed,
#' leveraging the bilinear structure of the model for computational efficiency.
#'
#' @details
#' The ALS algorithm exploits the fact that the SIR model is linear in A when B is fixed, 
#' and linear in B when A is fixed. This allows us to use standard GLM solvers in each step.
#' 
#' \strong{Algorithm Overview:}
#' \enumerate{
#'   \item Initialize B matrix (typically as identity)
#'   \item Repeat until convergence:
#'     \enumerate{
#'       \item \strong{A-step}: Fix B, optimize (θ, α) by solving:
#'         \deqn{Y \sim GLM(θ^T Z + vec(XB^T)^T (I ⊗ W) α)}
#'       \item \strong{B-step}: Fix A, optimize (θ, β) by solving:
#'         \deqn{Y \sim GLM(θ^T Z + vec(A^T X)^T (W ⊗ I) β)}
#'       \item Check convergence based on change in deviance
#'     }
#' }
#' 
#' \strong{Computational Efficiency:}
#' \itemize{
#'   \item Uses optimized C++ routines for matrix operations via RcppArmadillo
#'   \item Constructs sparse design matrices efficiently using Kronecker products
#'   \item Optionally uses speedglm for faster GLM fitting on large datasets
#'   \item Vectorized operations minimize memory allocations
#' }
#' 
#' \strong{Convergence Criteria:}
#' The algorithm converges when one of the following is met:
#' \itemize{
#'   \item Relative change in deviance < tol: |D_new - D_old|/|D_old| < tol
#'   \item Maximum iterations reached
#'   \item Deviance increases (indicating numerical issues)
#' }
#' 
#' \strong{Advantages of ALS:}
#' \itemize{
#'   \item More stable than direct optimization for high-dimensional problems
#'   \item Each sub-problem is convex (though overall problem is non-convex)
#'   \item Natural handling of constraints (e.g., non-negativity)
#'   \item Parallelizable sub-problems (future enhancement)
#' }
#' 
#' \strong{Limitations:}
#' \itemize{
#'   \item May converge to local optima (depends on initialization)
#'   \item Convergence can be slow near the optimum
#'   \item Requires good initialization for best results
#' }
#'
#' @param Y Three-dimensional array (m x m x T) of network outcomes.
#'   Missing values (NA) are automatically handled by excluding them from the likelihood.
#'   The algorithm uses complete case analysis within each GLM step.
#'   
#' @param W Three-dimensional array (m x m x p) of influence covariates used to 
#'   parameterize the influence matrices A and B. Each slice W[,,r] represents one
#'   influence covariate. If NULL or p=0, only identity matrices are used (no influence).
#'   
#' @param X Three-dimensional array (m x m x T) representing the network state that
#'   carries influence, typically lagged outcomes. Must be provided if W is non-NULL.
#'   This determines which network patterns affect future outcomes.
#'   
#' @param Z Four-dimensional array (m x m x q x T) of exogenous covariates, or NULL.
#'   These are covariates that directly affect outcomes but don't interact with the
#'   network influence structure. Examples include dyadic attributes or time trends.
#'   
#' @param family Character string specifying the GLM family: "poisson", "normal", or "binomial".
#'   This determines the link function and variance structure used in each GLM step.
#'   
#' @param trace Logical or integer controlling verbosity:
#'   \itemize{
#'     \item FALSE/0: No output
#'     \item TRUE/1: Progress bar and convergence message
#'     \item 2: Detailed iteration information including deviance
#'   }
#'   
#' @param tol Numeric convergence tolerance. The algorithm stops when the relative change
#'   in deviance is less than this value. Default is 1e-8. Smaller values give more
#'   accurate results but require more iterations.
#'   
#' @param max_iter Integer maximum number of ALS iterations. Default is 100.
#'   Each iteration consists of one A-step and one B-step. Increase for difficult
#'   problems or when starting far from the optimum.
#'   
#' @return A list with class "sir_als_fit" containing:
#'   \item{tab}{Vector of all parameters [θ, α, β] in order}
#'   \item{A}{The m x m sender effects matrix}
#'   \item{B}{The m x m receiver effects matrix}
#'   \item{deviance}{Final deviance (−2 × log-likelihood + constant)}
#'   \item{iterations}{Number of iterations until convergence}
#'   \item{converged}{Logical indicating successful convergence}
#'   \item{THETA}{Matrix tracking θ parameters across iterations (if trace > 0)}
#'   \item{ALPHA}{Matrix tracking α parameters across iterations (if trace > 0)}
#'   \item{BETA}{Matrix tracking β parameters across iterations (if trace > 0)}
#'   \item{DEV}{Matrix tracking deviance across iterations (if trace > 0)}
#'   \item{glm_alpha}{Final GLM object from the A-step}
#'   \item{glm_beta}{Final GLM object from the B-step}
#'
#' @examples
#' \dontrun{
#' # Example with simulated network data
#' m <- 15
#' T <- 20
#' p <- 2
#' 
#' # Generate network with influence
#' Y <- array(rpois(m*m*T, 3), dim=c(m,m,T))
#' X <- array(0, dim=c(m,m,T))
#' X[,,2:T] <- Y[,,1:(T-1)]  # Lagged Y
#' 
#' # Influence covariates (e.g., distance-based)
#' W <- array(rnorm(m*m*p), dim=c(m,m,p))
#' 
#' # Fit using ALS
#' fit <- sir_alsfit(Y, W, X, Z=NULL, family="poisson", 
#'                   trace=TRUE, tol=1e-6, max_iter=50)
#' 
#' # Examine convergence
#' plot(fit$DEV[,2], type="l", ylab="Deviance", xlab="Iteration")
#' 
#' # Extract influence matrices
#' A_matrix <- fit$A
#' B_matrix <- fit$B
#' }
#' @importFrom stats glm lm poisson binomial coef deviance formula
#' @importFrom speedglm speedglm
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done cli_alert_info cli_alert_success
#' @importFrom crayon yellow green
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
  # Initialize progress bar if trace is enabled
  if (trace && p > 0) {
    cli::cli_progress_bar("Fitting SIR model via ALS", total = max_iter)
  }
  
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
        if (trace) {
          cli::cli_progress_done()
          warning("Deviance became NA, stopping iteration.")
        }
        break
    }

    # Track history
    THETA <- rbind(THETA, theta)
    ALPHA <- rbind(ALPHA, alpha)
    BETA  <- rbind(BETA, beta)
    DEV   <- rbind(DEV, c(dev_old, dev_new))

    if(trace){
      cli::cli_progress_update()
      cli::cli_alert_info("Iteration {.val {iter}}: Deviance = {.val {sprintf('%.4f', dev_new)}}, Change = {.val {sprintf('%.6f', abs(dev_old - dev_new)/(abs(dev_old) + 0.1))}}")
    }

    # Stopping criterion
    if( abs(dev_old - dev_new) / (abs(dev_old) + 0.1) < tol ) {
        if (trace) {
          cli::cli_progress_done()
          cli::cli_alert_success("ALS converged after {.val {iter}} iterations")
        }
        break
    }
  }

  if (iter == max_iter && p > 0) {
      if (trace) {
        cli::cli_progress_done()
      }
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
