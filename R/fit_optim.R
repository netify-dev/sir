
#' Fit SIR Model via Direct Optimization
#'
#' @description
#' Implements direct optimization of the Social Influence Regression model using the 
#' Broyden-Fletcher-Goldfarb-Shanno (BFGS) quasi-Newton method. This approach optimizes
#' all parameters simultaneously, treating the full likelihood as a single objective function.
#'
#' @details
#' Unlike the ALS method which alternates between parameter blocks, direct optimization
#' treats the SIR model as a single non-linear optimization problem:
#' 
#' \deqn{\min_{\theta, \alpha, \beta} -\ell(\theta, \alpha, \beta | Y, X, W, Z)}
#' 
#' where \eqn{\ell} is the log-likelihood function for the specified family.
#' 
#' \strong{Optimization Strategy:}
#' \itemize{
#'   \item Uses BFGS algorithm with analytical gradients computed via C++
#'   \item Gradients leverage the chain rule through the bilinear structure
#'   \item Line search ensures sufficient decrease in objective
#'   \item Hessian approximation updated using gradient information
#' }
#' 
#' \strong{Gradient Computation:}
#' The gradient with respect to parameters is:
#' \itemize{
#'   \item \eqn{\nabla_\theta \ell}: Direct gradient for exogenous covariates
#'   \item \eqn{\nabla_\alpha \ell}: Gradient flows through A = alpha_0 I + sum_r alpha_r W_r
#'   \item \eqn{\nabla_\beta \ell}: Gradient flows through B = beta_0 I + sum_r beta_r W_r
#' }
#' 
#' \strong{Initialization:}
#' If no starting values provided, uses smart initialization:
#' \enumerate{
#'   \item Fit GLM ignoring network effects to initialize theta
#'   \item Set alpha and beta to small random values or zeros
#'   \item Scale parameters based on outcome variance
#' }
#' 
#' \strong{Advantages over ALS:}
#' \itemize{
#'   \item Can find better global optima (not restricted to alternating updates)
#'   \item Faster convergence when near optimum
#'   \item Provides Hessian for standard error calculation
#'   \item Single convergence criterion
#' }
#' 
#' \strong{Disadvantages:}
#' \itemize{
#'   \item Less stable for high-dimensional problems (p large)
#'   \item Can fail with poor initialization
#'   \item Memory intensive for large networks
#'   \item Sensitive to scaling of covariates
#' }
#' 
#' \strong{Convergence Diagnostics:}
#' The optimizer reports convergence codes:
#' \itemize{
#'   \item 0: Successful convergence
#'   \item 1: Maximum iterations reached
#'   \item 10: Degeneracy of Nelder-Mead simplex
#'   \item 51: Warning from L-BFGS-B
#'   \item 52: Error from L-BFGS-B
#' }
#'
#' @param Y Three-dimensional array (m x m x T) of network outcomes.
#'   Missing values (NA) are handled by excluding them from the likelihood.
#'   Large networks (m > 100) may cause memory issues.
#'   
#' @param W Three-dimensional array (m x m x p) of influence covariates.
#'   Each slice parameterizes the influence matrices. If NULL, no network
#'   influence structure is included (reduces to standard GLM).
#'   
#' @param X Three-dimensional array (m x m x T) representing the network state
#'   carrying influence. Typically lagged Y. Required if W is provided.
#'   
#' @param Z Four-dimensional array (m x m x q x T) of exogenous covariates.
#'   These enter the model linearly without network interactions.
#'   
#' @param family Character string: "poisson", "normal", or "binomial".
#'   Determines the likelihood and link function used.
#'   
#' @param trace Integer controlling optimizer output:
#'   \itemize{
#'     \item 0: No output (default)
#'     \item 1: Final convergence report
#'     \item 2: Progress at each iteration
#'     \item 3+: Detailed debugging information
#'   }
#'   
#' @param start Optional numeric vector of starting values [θ, α, β].
#'   Length must equal q + 2p. If NULL, uses smart initialization.
#'   Good starting values dramatically improve convergence.
#'   
#' @return A list with class "sir_optim_fit" containing:
#'   \item{tab}{Vector of optimized parameters [θ, α, β]}
#'   \item{A}{The m x m sender effects matrix}
#'   \item{B}{The m x m receiver effects matrix}
#'   \item{convergence}{Convergence code from optim (0 = success)}
#'   \item{message}{Convergence message from optimizer}
#'   \item{iterations}{Number of function evaluations}
#'   \item{value}{Final negative log-likelihood}
#'   \item{hessian}{Approximate Hessian at optimum (if requested)}
#'   \item{gradient}{Final gradient (should be near zero)}
#'   
#' @examples
#' \dontrun{
#' # Small network example
#' m <- 10
#' T <- 15
#' p <- 2
#' 
#' Y <- array(rpois(m*m*T, 2), dim=c(m,m,T))
#' X <- array(0, dim=c(m,m,T))
#' X[,,2:T] <- Y[,,1:(T-1)]
#' W <- array(rnorm(m*m*p, sd=0.1), dim=c(m,m,p))
#' 
#' # Fit with direct optimization
#' fit_optim <- sir_optfit(Y, W, X, Z=NULL,
#'                         family="poisson", trace=1)
#'                         
#' # Check convergence
#' if(fit_optim$convergence == 0) {
#'   cat("Optimization successful\n")
#'   print(fit_optim$tab)
#' }
#' 
#' # Compare with ALS
#' fit_als <- sir_alsfit(Y, W, X, Z=NULL,
#'                       family="poisson")
#' 
#' # Often similar results but different paths
#' }
#' @importFrom stats optim lm glm poisson binomial coef formula
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_warning
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
  if (trace > 0) {
    cli::cli_alert_info("Starting BFGS optimization with {.val {length(start)}} parameters")
  }
  
  fit <- optim(par=start, fn=objfun, gr=gradfun, method="BFGS", control=list(trace=trace, maxit=500))
  
  if (trace > 0) {
    if (fit$convergence == 0) {
      cli::cli_alert_success("Optimization converged: NLL = {.val {sprintf('%.4f', fit$value)}}, {.val {fit$counts[1]}} iterations")
    } else {
      cli::cli_alert_warning("Optimization did not converge (code {.val {fit$convergence}})")
    }
  }

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
