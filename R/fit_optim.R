
#' Fit SIR model via Direct Optimization (BFGS)
#'
#' Minimizes the Negative Log-Likelihood (NLL) using BFGS, utilizing the analytical gradient calculated via Rcpp.
#'
#' @param Y (m x m x T) outcome array.
#' @param W (m x m x p) influence covariates.
#' @param X (m x m x T) bilinear covariates.
#' @param Z (m x m x q x T) exogenous covariates.
#' @param family Distribution family.
#' @param trace Integer level for optim tracing.
#' @param start Optional starting parameter vector.
#' @return A list containing the estimated parameters and optimization results.
#' @importFrom stats optim lm glm poisson binomial coef formula
#' @export
sir_optfit <- function(Y, W, X, Z, family, trace=0, start=NULL) {
  p <- if (is.null(W)) 0 else dim(W)[3]
  q <- if (is.null(Z)) 0 else dim(Z)[3]

  # ---- Initialization ----
  if(is.null(start)) {
    # Naive start: GLM ignoring bilinear portion
    Y_flat <- flatten_Y(Y)
    Z_flat <- flatten_Z(Z)

    if (q > 0) {
        glmData <- data.frame(Y = Y_flat, Z_flat)
        form <- formula(paste0("Y ~ -1 + ", paste(colnames(Z_flat), collapse=" + ")))

        if (family == "normal") {
            fit0 <- lm(form, data=glmData)
        } else {
            family_obj <- switch(family,
                                "poisson" = poisson(),
                                "binomial" = binomial())
            fit0 <- glm(form, data=glmData, family=family_obj)
        }
        theta <- coef(fit0)
        # Handle potential NA coefficients
        theta[is.na(theta)] <- 0
    } else {
        theta <- numeric(0)
    }

    # Small random for alpha[-1], beta
    if (p > 0) {
        # Initialize with slightly larger values for stability in optim
        start <- c(theta, rnorm(p-1, sd=0.05), rnorm(p, sd=0.05))
    } else {
        start <- theta
    }
  }

  # Prepare Z list for C++ (Handles 3D/4D Z correctly)
  Z_list <- prepare_Z_list(Z)

  # Define objective function (NLL)
  objfun <- function(par){
    # Use the R implementation for NLL calculation (fast enough and uses R stats functions)
    # We pass Z (the original R array) to mll_sir/eta_tab
    mll_sir(par, Y, W, X, Z, family)
  }

  # Define gradient function (using Rcpp backend)
  gradfun <- function(par){
    # Use the C++ implementation for the gradient of NLL
    # We pass Z_list (the list of cubes) to cpp_mll_gH
    gH <- cpp_mll_gH(par, Y, W, X, Z_list, family)
    as.numeric(gH$grad)
  }

  # Run optimization
  # We minimize the NLL using the gradient of the NLL.
  fit <- optim(par=start, fn=objfun, gr=gradfun, method="BFGS", control=list(trace=trace, maxit=500))

  tab <- fit$par

  # Parse final parameters
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
      } else {
          a <- numeric(0)
      }
      b_start = q + p
      b <- tab[b_start:length(tab)]
  } else {
      a <- numeric(0)
      b <- numeric(0)
  }

  list(
    theta=theta,
    a=a,
    b=b,
    tab=tab,
    value=fit$value,
    convergence=fit$convergence,
    counts=fit$counts,
    iterations=fit$counts[1] # Use function evaluations as iteration count proxy
  )
}
