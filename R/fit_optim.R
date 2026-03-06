
#' Fit SIR Model via Direct Optimization
#'
#' @description
#' Optimizes all SIR model parameters simultaneously using BFGS with
#' analytical gradients computed via C++. An alternative to the ALS method
#' that can be faster for small problems or when good starting values are
#' available.
#'
#' @details
#' This function treats the SIR negative log-likelihood as a single
#' optimization problem over all parameters (theta, alpha, beta). It uses
#' the BFGS quasi-Newton method with analytical gradients from the C++
#' backend.
#'
#' If no starting values are provided, theta is initialized by fitting a GLM
#' ignoring network effects, and alpha/beta are set to small random values.
#'
#' Convergence code 0 indicates success; code 1 means the maximum number of
#' iterations was reached. Only supports static (3D) W.
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
#' @param start Optional numeric vector of starting values [theta, alpha, beta].
#'   Length must equal q + 2p. If NULL, uses smart initialization.
#'   Good starting values dramatically improve convergence.
#'   
#' @return A list with class "sir_optim_fit" containing:
#'   \item{tab}{Vector of optimized parameters [theta, alpha, beta]}
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

  # initialization: GLM ignoring bilinear portion
  if(is.null(start)) {
	Y_flat <- flatten_Y(Y)
	Z_flat <- flatten_Z(Z)

	if (q > 0) {
		glm_data <- data.frame(Y = Y_flat, Z_flat)
		form <- formula(paste0("Y ~ -1 + ", paste(colnames(Z_flat), collapse=" + ")))

		if (family == "normal") {
			fit0 <- lm(form, data=glm_data)
		} else {
			family_obj <- switch(family,
								"poisson" = poisson(),
								"binomial" = binomial())
			fit0 <- glm(form, data=glm_data, family=family_obj)
		}
		theta <- coef(fit0)
		# handle NA coefficients
		theta[is.na(theta)] <- 0
	} else {
		theta <- numeric(0)
	}

	# small random start for alpha[-1], beta
	if (p > 0) {
		start <- c(theta, rnorm(p-1, sd=0.05), rnorm(p, sd=0.05))
	} else {
		start <- theta
	}
  }

  # prepare Z list for C++
  Z_list <- prepare_Z_list(Z)

  # objective: negative log-likelihood
  objfun <- function(par){
	mll_sir(par, Y, W, X, Z, family)
  }

  # gradient via C++ backend
  gradfun <- function(par){
	gH <- cpp_mll_gH(par, Y, W, X, Z_list, family)
	as.numeric(gH$grad)
  }

  # run BFGS optimization
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

  # parse final parameters
  if (q > 0) {
	  theta <- tab[1:q]
  } else {
	  theta <- numeric(0)
  }

  if (p > 0) {
	  if (p > 1) {
		  a_start <- q + 1
		  a_end <- q + p - 1
		  a <- tab[a_start:a_end]
	  } else {
		  a <- numeric(0)
	  }
	  b_start <- q + p
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
	iterations=fit$counts[1]
  )
}
