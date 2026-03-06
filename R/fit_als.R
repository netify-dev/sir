
#' Fit SIR Model via Alternating Least Squares (ALS)
#'
#' @description
#' Fits the SIR model by alternating between optimizing sender effects (alpha)
#' with receiver effects fixed, and vice versa. Each sub-step is a standard
#' GLM, making this approach more stable than direct optimization for
#' high-dimensional problems.
#'
#' @details
#' The algorithm exploits the bilinear structure: when B is fixed, the model
#' is linear in alpha (and theta), and vice versa. Each iteration solves two
#' GLMs using \code{speedglm} (if available) or \code{stats::glm}.
#'
#' Convergence is declared when the relative change in deviance falls below
#' \code{tol}, or when \code{max_iter} is reached. Supports both static (3D)
#' and dynamic (4D) W arrays.
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
#' @param fix_receiver Logical. If TRUE, fixes B = I (identity matrix) and
#'   estimates only (theta, alpha) via a single GLM step. This eliminates the
#'   bilinear identification problem by removing the receiver influence channel.
#'   The model reduces to a standard GLM with well-conditioned standard errors.
#'   Default is FALSE.
#'
#' @param kron_mode Logical. If TRUE, replaces separate (alpha, beta) with a
#'   single p x p coefficient matrix C. Not yet implemented. Default is FALSE.
#'
#' @param dynamic_W Logical. If TRUE, W is treated as a 4D array
#'   (m x m x p x T) with time-varying influence covariates. Design matrices
#'   are constructed in R rather than C++. Default is FALSE.
#'
#' @return A list with class "sir_als_fit" containing:
#'   \item{tab}{Vector of all parameters [theta, alpha, beta] in order}
#'   \item{A}{The m x m sender effects matrix}
#'   \item{B}{The m x m receiver effects matrix}
#'   \item{deviance}{Final deviance (-2 * log-likelihood + constant)}
#'   \item{iterations}{Number of iterations until convergence}
#'   \item{converged}{Logical indicating successful convergence}
#'   \item{THETA}{Matrix tracking theta parameters across iterations}
#'   \item{ALPHA}{Matrix tracking alpha parameters across iterations}
#'   \item{BETA}{Matrix tracking beta parameters across iterations}
#'   \item{DEV}{Matrix tracking deviance across iterations}
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
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done cli_alert_info cli_alert_success
#' @export
sir_alsfit <- function(Y, W, X, Z, family, trace=FALSE, tol=1e-8, max_iter=100,
					   fix_receiver=FALSE, kron_mode=FALSE, dynamic_W=FALSE) {
  p <- if (is.null(W)) 0 else dim(W)[3]
  q <- if (is.null(Z)) 0 else dim(Z)[3]
  n1 <- dim(Y)[1]
  n2 <- dim(Y)[2]
  m <- n1  # backward compat for square case
  T_len <- dim(Y)[3]
  N_flat <- n1 * n2 * T_len  # total number of entries in flattened Y

  # pre-convert 4D W to list-of-cubes for C++ when dynamic
  W_field <- if (dynamic_W && p > 0) prepare_W_field(W) else NULL

  # determine GLM function and family object
  if (family == "normal") {
	glm_fun <- stats::lm
	family_obj <- NULL
	use_speedglm <- FALSE
  } else {
	# use speedglm if available
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
						 cli::cli_abort("Unsupported family for GLM: {.val {family}}."))
  }

  # initialization
  Y_flat <- flatten_Y(Y)
  Z_flat <- flatten_Z(Z)

  if (q > 0) {
	glm_data <- data.frame(Y = Y_flat, Z_flat)
	form <- formula(paste0("Y ~ -1 + ", paste(colnames(Z_flat), collapse=" + ")))

	if (family == "normal") {
		fit0 <- glm_fun(form, data=glm_data)
	} else {
		fit0 <- glm_fun(form, data=glm_data, family=family_obj)
	}
	theta <- coef(fit0)
	# handle NA coefficients from initialization
	theta[is.na(theta)] <- 0
	dev_new <- if (family == "normal") sum(fit0$residuals^2) else deviance(fit0)

  } else {
	  theta <- numeric(0)
	  # null deviance when q=0
	  mean_Y <- mean(Y_flat, na.rm=TRUE)
	  if (family == "normal") {
		  dev_new <- sum((Y_flat - mean_Y)^2, na.rm=TRUE)
	  } else if (family == "poisson") {
		  # avoid log(0)
		  if (mean_Y <= 0) mean_Y <- 1e-6
		  dev_new <- 2 * sum(ifelse(Y_flat > 0, Y_flat * log(Y_flat/mean_Y), 0) - (Y_flat - mean_Y), na.rm=TRUE)
	  } else if (family == "binomial") {
		  p0 <- mean_Y
		  if (p0 <= 0) p0 <- 1e-6
		  if (p0 >= 1) p0 <- 1 - 1e-6
		  dev_new <- -2 * sum(Y_flat * log(p0) + (1-Y_flat) * log(1-p0), na.rm=TRUE)
	  }
  }


  # initialize alpha, beta with small random values
  if (p > 0) {
	  alpha <- rnorm(p, sd=0.05)
	  beta  <- rnorm(p, sd=0.05)
  } else {
	  alpha <- numeric(0)
	  beta <- numeric(0)
  }


  # track iteration history
  dev_old <- Inf
  THETA <- matrix(theta, nrow=1)
  ALPHA <- matrix(alpha, nrow=1)
  BETA  <- matrix(beta, nrow=1)
  DEV   <- matrix(c(dev_old, dev_new), nrow=1)

  # fix_receiver: single GLM with B = identity
  if (fix_receiver && p > 0) {
	if (kron_mode) cli::cli_abort("Cannot use both {.arg fix_receiver} and {.arg kron_mode}.")

	# build design matrix: column k = vec(W_k %*% X_t) across all t
	Wfix_flat <- matrix(0, N_flat, p)
	for (k in 1:p) {
	  wk_x <- array(0, dim = c(n1, n2, T_len))
	  for (t in 1:T_len) {
		W_kt <- if (dynamic_W) W[,,k,t] else W[,,k]
		wk_x[,,t] <- W_kt %*% X[,,t]
	  }
	  Wfix_flat[, k] <- c(wk_x)
	}
	colnames(Wfix_flat) <- paste0("WA", 1:p)

	if (q > 0) {
	  X_design <- cbind(Z_flat, Wfix_flat)
	} else {
	  X_design <- Wfix_flat
	}

	glm_data <- data.frame(Y = Y_flat, X_design)
	form <- formula(paste0("Y ~ -1 + ", paste(colnames(X_design), collapse = " + ")))

	if (family == "normal") {
	  fit_fix <- glm_fun(form, data = glm_data)
	} else {
	  fit_fix <- glm_fun(form, data = glm_data, family = family_obj)
	}

	co <- coef(fit_fix)
	co[is.na(co)] <- 0

	if (q > 0) {
	  theta <- co[1:q]
	  alpha <- co[-(1:q)]
	} else {
	  theta <- numeric(0)
	  alpha <- co
	}

	dev_val <- if (family == "normal") sum(fit_fix$residuals^2) else deviance(fit_fix)

	if (trace) {
	  cli::cli_alert_success("fix_receiver: single GLM converged (deviance = {.val {sprintf('%.4f', dev_val)}})")
	}

	return(list(
	  theta = theta,
	  a = alpha,
	  b = numeric(0),
	  tab = c(theta, alpha),
	  ALPHA = matrix(alpha, nrow = 1),
	  BETA = matrix(nrow = 0, ncol = 0),
	  THETA = matrix(theta, nrow = 1),
	  DEV = matrix(c(Inf, dev_val), nrow = 1),
	  iterations = 1,
	  fix_receiver = TRUE,
	  glm_fit = fit_fix
	))
  }

  # kron_mode placeholder
  if (kron_mode) cli::cli_abort("{.arg kron_mode} is not yet implemented.")

  # iterative ALS updates
  iter <- 0L
  if (trace && p > 0) {
	cli::cli_progress_bar("Fitting SIR model via ALS", total = max_iter)
  }

  for (iter in 1:max_iter) {

	if (p == 0) break

	# fix beta, update (theta, alpha)
	# construct Wbeta design matrix
	if (dynamic_W) {
		Wbeta_flat <- cpp_construct_Wbeta_design_dyn(W_field, X, beta)
	} else {
		Wbeta_flat <- cpp_construct_Wbeta_design(W, X, beta)
	}
	colnames(Wbeta_flat) <- paste0("WB", 1:p)

	if (q > 0) {
		X_design <- cbind(Z_flat, Wbeta_flat)
	} else {
		X_design <- Wbeta_flat
	}
	glm_data <- data.frame(Y = Y_flat, X_design)
	form1 <- formula(paste0("Y ~ -1 + ", paste(colnames(X_design), collapse=" + ")))

	start_vals <- c(theta, alpha)
	fit_a <- fit_glm_wrapper(form1, glm_data, family, family_obj, glm_fun, use_speedglm, start_vals, trace)

	co_a  <- coef(fit_a)
	co_a[is.na(co_a)] <- 0

	if (q > 0) {
		theta <- co_a[1:q]
		alpha <- co_a[-(1:q)]
	} else {
		alpha <- co_a
	}

	# fix alpha, update (theta, beta)
	# construct Walpha design matrix
	if (dynamic_W) {
		Walpha_flat <- cpp_construct_Walpha_design_dyn(W_field, X, alpha)
	} else {
		Walpha_flat <- cpp_construct_Walpha_design(W, X, alpha)
	}
	colnames(Walpha_flat) <- paste0("WA", 1:p)

	if (q > 0) {
		X_design <- cbind(Z_flat, Walpha_flat)
	} else {
		X_design <- Walpha_flat
	}
	glm_data <- data.frame(Y = Y_flat, X_design)
	form2 <- formula(paste0("Y ~ -1 + ", paste(colnames(X_design), collapse=" + ")))

	start_vals <- c(theta, beta)
	fit_b <- fit_glm_wrapper(form2, glm_data, family, family_obj, glm_fun, use_speedglm, start_vals, trace)

	co_b  <- coef(fit_b)
	co_b[is.na(co_b)] <- 0

	if (q > 0) {
		theta <- co_b[1:q]
		beta  <- co_b[-(1:q)]
	} else {
		beta <- co_b
	}

	# check deviance
	dev_old <- dev_new
	dev_new <- if (family == "normal") sum(fit_b$residuals^2) else deviance(fit_b)

	# handle NA deviance
	if (is.na(dev_new)) {
		if (trace) {
		  cli::cli_progress_done()
		}
		cli::cli_warn("Deviance became NA at iteration {.val {iter}}, stopping.")
		dev_new <- dev_old
		break
	}

	# track history
	THETA <- rbind(THETA, theta)
	ALPHA <- rbind(ALPHA, alpha)
	BETA  <- rbind(BETA, beta)
	DEV   <- rbind(DEV, c(dev_old, dev_new))

	if(trace){
	  cli::cli_progress_update()
	  cli::cli_alert_info("Iteration {.val {iter}}: Deviance = {.val {sprintf('%.4f', dev_new)}}, Change = {.val {sprintf('%.6f', abs(dev_old - dev_new)/(abs(dev_old) + 0.1))}}")
	}

	# stopping criterion
	if( abs(dev_old - dev_new) / (abs(dev_old) + 0.1) < tol ) {
		if (trace) {
		  cli::cli_progress_done()
		  cli::cli_alert_success("ALS converged after {.val {iter}} iterations")
		}
		break
	}
  }

  # determine convergence
  als_converged <- (p == 0) ||
	  (abs(dev_old - dev_new) / (abs(dev_old) + 0.1) < tol)

  if (!als_converged && p > 0) {
	  if (trace) {
		cli::cli_progress_done()
	  }
	  cli::cli_warn("ALS did not converge within the maximum number of iterations.")
  }

  # impose alpha_1=1 constraint and rescale
  if (p > 0) {
	  if (abs(alpha[1]) < 1e-8) {
		  cli::cli_warn("Estimated alpha[1] is near zero. Rescaling might be unstable. Using pseudoinverse stabilization.")
		  # stabilize near-zero alpha[1]
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
	iterations = iter,
	converged = als_converged
  )
}

# helper for GLM fitting with fallback on failure
fit_glm_wrapper <- function(form, data, family, family_obj, glm_fun, use_speedglm, start_vals, trace) {
	if (family == "normal") {
		fit <- glm_fun(form, data=data)
	} else {
		# ensure start values match predictor count
		if (length(start_vals) != ncol(data)-1) {
			 if (trace) cli::cli_warn("Start values length mismatch. Using default initialization.")
			 fit <- glm_fun(form, data=data, family=family_obj)
		} else {
			# try start values, fall back on failure
			tryCatch({
				fit <- glm_fun(form, data=data, family=family_obj, start=start_vals)
			}, error = function(e) {
				if (trace) cli::cli_warn("GLM fitting failed with start values ({e$message}), falling back to default initialization.")
				# fallback without start values
				tryCatch({
					 fit <<- glm_fun(form, data=data, family=family_obj)
				}, error = function(e2) {
					# final fallback to stats::glm
					if (use_speedglm) {
						if (trace) cli::cli_warn("{.pkg speedglm} failed entirely, falling back to {.fn stats::glm}.")
						fit <<- stats::glm(form, data=data, family=family_obj)
					} else {
						cli::cli_abort("{.fn stats::glm} failed during ALS iteration: {e2$message}")
					}
				})
			})
		}
	}
	return(fit)
}
