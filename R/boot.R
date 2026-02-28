
#' Bootstrap Standard Errors for SIR Model
#'
#' @description
#' Computes bootstrap standard errors and confidence intervals for SIR model
#' parameters. Uses block bootstrap (resampling time periods) or parametric
#' bootstrap (simulating from the fitted model) to avoid reliance on the
#' Hessian matrix, which can be singular or ill-conditioned.
#'
#' @param sir_fit A fitted sir object from \code{\link{sir}}.
#' @param R Integer. Number of bootstrap replicates (default 200).
#' @param type Character. Bootstrap type: "block" resamples time periods,
#'   "parametric" simulates new data from the fitted model.
#' @param seed Optional integer random seed for reproducibility.
#' @param trace Logical. Print progress every 10 replicates.
#' @return A list containing:
#' \describe{
#'   \item{coefs}{R x n_params matrix of bootstrap coefficient estimates}
#'   \item{se}{Bootstrap standard errors (column SDs)}
#'   \item{ci_lo}{Lower 2.5\% percentile confidence bounds}
#'   \item{ci_hi}{Upper 97.5\% percentile confidence bounds}
#'   \item{n_valid}{Number of successful bootstrap replicates}
#'   \item{n_total}{Total number attempted}
#' }
#' @export
boot_sir <- function(sir_fit, R = 200, type = c("block", "parametric"),
                     seed = NULL, trace = FALSE) {
    type <- match.arg(type)
    if (!is.null(seed)) set.seed(seed)

    Y <- sir_fit$Y
    W <- sir_fit$W
    X <- sir_fit$X
    Z <- sir_fit$Z
    family <- sir_fit$family
    T_len <- dim(Y)[3]
    n_params <- length(sir_fit$tab)

    boot_coefs <- matrix(NA, R, n_params)
    colnames(boot_coefs) <- rownames(sir_fit$summ)

    for (b in 1:R) {
        if (trace && b %% 10 == 0) cat(sprintf("Bootstrap %d/%d\n", b, R))

        if (type == "block") {
            # Resample time periods with replacement
            t_idx <- sample(1:T_len, T_len, replace = TRUE)
            Y_b <- Y[,, t_idx, drop = FALSE]
            X_b <- X[,, t_idx, drop = FALSE]
            Z_b <- if (!is.null(Z) && length(dim(Z)) == 4) {
                Z[,,, t_idx, drop = FALSE]
            } else {
                Z
            }
        } else {
            # Parametric: simulate from fitted model
            eta <- eta_tab(sir_fit$tab, W, X, Z)
            if (family == "normal") {
                sigma <- sqrt(sir_fit$sigma2)
                Y_b <- array(rnorm(length(eta), mean = eta, sd = sigma),
                             dim = dim(Y))
            } else if (family == "poisson") {
                lambda <- exp(eta)
                lambda[lambda > 1e6] <- 1e6
                Y_b <- array(rpois(length(eta), lambda = lambda),
                             dim = dim(Y))
            } else {
                prob <- 1 / (1 + exp(-eta))
                Y_b <- array(rbinom(length(eta), 1, prob),
                             dim = dim(Y))
            }
            # Preserve diagonal NA structure
            for (tt in 1:T_len) diag(Y_b[,,tt]) <- NA
            X_b <- X
            Z_b <- Z
        }

        tryCatch({
            fit_b <- sir(Y_b, W, X_b, Z_b, family = family,
                         method = "ALS", calcSE = FALSE,
                         fix_receiver = isTRUE(sir_fit$fix_receiver),
                         kron_mode = isTRUE(sir_fit$kron_mode),
                         max_iter = 100, tol = 1e-6)
            boot_coefs[b, ] <- fit_b$tab
        }, error = function(e) {
            # Leave as NA for failed replicates
        })
    }

    # Compute statistics from successful replicates
    valid <- apply(boot_coefs, 1, function(r) !any(is.na(r)))
    n_valid <- sum(valid)

    if (n_valid < 10) {
        warning(sprintf(
            "Only %d valid bootstrap replicates (of %d). Results unreliable.",
            n_valid, R))
    }

    boot_se <- apply(boot_coefs[valid, , drop = FALSE], 2, sd)
    boot_ci <- apply(boot_coefs[valid, , drop = FALSE], 2,
                     quantile, probs = c(0.025, 0.975))

    list(
        coefs = boot_coefs,
        se = boot_se,
        ci_lo = boot_ci[1, ],
        ci_hi = boot_ci[2, ],
        n_valid = n_valid,
        n_total = R
    )
}
