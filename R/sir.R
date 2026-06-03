
#' Social Influence Regression (SIR) Model
#'
#' @description
#' Fits a Social Influence Regression model for network data with lagged,
#' model-implied association channels. The SIR model captures how prior network
#' states predict later outcomes through bilinear interaction terms, allowing for
#' both sender and receiver channels in directed networks. Causal interpretation
#' requires additional research-design assumptions.
#' 
#' The model decomposes network influence into two components:
#' \itemize{
#'   \item \strong{Sender influence (A matrix)}: A[i,k] measures how much node k's
#'     behavior (via X) shapes node i's outgoing ties.
#'   \item \strong{Receiver influence (B matrix)}: B[j,l] measures how node l's
#'     position shapes node j's incoming ties.
#' }
#'
#' @details
#' The SIR model specifies the expected outcome for the directed edge from node i to node j 
#' at time t as:
#' 
#' \deqn{g(\mu_{i,j,t}) = \theta^T z_{i,j,t} + \sum_{k,l} X_{k,l,t} A_{i,k} B_{j,l}}
#'
#' Where:
#' \itemize{
#'   \item \eqn{g(\cdot)} is the family link function (identity for normal, log
#'     for poisson, logit for binomial). The exogenous and bilinear terms enter
#'     the \emph{linear predictor}, not the mean directly, so for poisson and
#'     binomial they act on the log / logit scale.
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
#'   \item \eqn{A = \sum_{r=1}^{p} \alpha_r W_r} (sender effects, \eqn{\alpha_1 = 1} fixed)
#'   \item \eqn{B = \sum_{r=1}^{p} \beta_r W_r} (receiver effects)
#' }
#' 
#' This parameterization reduces the number of parameters from \eqn{O(m^2)} to \eqn{O(p)}, where \eqn{p \ll m}.
#' 
#' @section Estimation Methods:
#' 
#' \strong{Alternating GLM/IRLS updates (method = "ALS"):}
#' \itemize{
#'   \item Iteratively updates the sender and receiver influence coefficients
#'     with GLM/IRLS subproblems
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
#' Both engines report the same public parameterization with
#' \eqn{\alpha_1 = 1}. BFGS optimizes that reduced vector directly; the default
#' alternating engine updates working alpha/beta coordinates with GLM/IRLS
#' subproblems and then normalizes to
#' \eqn{\alpha_1 = 1} when possible. The full bilinear log-likelihood is often
#' flat along a ridge: the alternating engine and BFGS typically reach the same
#' log-likelihood (to well within 1\%) yet can return parameter vectors that
#' differ by a few tenths on weakly identified influence coefficients. Compare fits by log-likelihood /
#' fitted values rather than by raw coefficients, and use
#' \code{fix_receiver = TRUE} to remove the alpha/beta scaling ambiguity when a
#' simpler sender-side channel is substantively adequate.
#'
#' @section Distribution Families:
#' 
#' \strong{Poisson:} For count data (e.g., number of interactions)
#' \itemize{
#'   \item Link function: log
#'   \item Variance function: \eqn{V(\mu) = \mu}
#'   \item Use when: Y_ijt represents counts
#' }
#'
#' \strong{Normal:} For continuous data (e.g., trade volumes, distances)
#' \itemize{
#'   \item Link function: identity
#'   \item Variance function: \eqn{V(\mu) = \sigma^2}
#'   \item Use when: Y_ijt is continuous and approximately normal
#' }
#'
#' \strong{Binomial:} For binary data (e.g., presence/absence of ties)
#' \itemize{
#'   \item Link function: logit
#'   \item Variance function: \eqn{V(\mu) = \mu(1 - \mu)}
#'   \item Use when: Y_ijt is binary (0/1)
#' }
#'
#' @param Y A three-dimensional array of dimensions (m x m x T) containing the network 
#'   outcomes. Y[i,j,t] represents the directed outcome from node i to node j at time t.
#'   Can contain NA values for missing observations. For one-mode square networks,
#'   self-tie diagonals are excluded from fitting, likelihood evaluation, and
#'   reported observation counts.
#'   
#' @param W Optional influence covariate array, either:
#'   \itemize{
#'     \item \strong{3D array} (m x m x p): Static influence covariates.
#'       W[i,k,r] represents the r-th covariate for the channel by which source
#'       node k can shape target node i; the same orientation is used for the
#'       receiver-side matrix B[j,l]. Directed W slices therefore need the same
#'       target-by-source orientation as A and B, not necessarily the same
#'       sender-to-receiver orientation as Y[i,j,t]. The same W is used for all
#'       time periods.
#'     \item \strong{4D array} (m x m x p x T): Dynamic (time-varying) influence
#'       covariates. W[i,k,r,t] allows the influence structure to change over time.
#'       Parameters (alpha, beta) are still estimated jointly across all periods,
#'       but the influence matrices A_t and B_t vary with t. Only ALS method is
#'       supported for 4D W.
#'   }
#'   Common choices include graph Laplacians, geographic distance matrices, or
#'   node-level covariates expanded to edge-level.
#'   For one-mode square networks, W diagonals are set to zero so self-influence
#'   channels are not estimated. If NULL or p=0, no network influence structure
#'   is included.
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
#' @param method Character string specifying the estimation method. Either
#'   \code{"ALS"} (the default alternating GLM/IRLS engine; the name is retained
#'   for API compatibility) or \code{"optim"} (direct optimization via BFGS).
#'   
#' @param calc_se Logical indicating whether to calculate standard errors for the
#'   parameters. Standard errors are computed using the observed information matrix.
#'   Setting to FALSE speeds up computation when uncertainty quantification is not needed.
#'
#' @param fix_receiver Logical. If TRUE, fixes B = I (identity matrix) and
#'   estimates only (theta, alpha). This eliminates the bilinear identification
#'   problem (scaling ambiguity between A and B) by removing the receiver
#'   influence channel. The model becomes a standard GLM, yielding model-based
#'   standard errors under the usual GLM assumptions. Appropriate when receiver
#'   effects are negligible.
#'   Default is FALSE.
#'
#' @param symmetric Logical. If TRUE, treats the network as undirected
#'   (symmetric). The function uses only upper-triangle observations for
#'   fitting and sets \code{fix_receiver = TRUE}, giving an upper-triangle,
#'   sender-side representation of undirected data rather than a fully
#'   order-invariant undirected bilinear model. Continuous asymmetric outcomes
#'   are averaged across upper and lower triangles; Poisson and Binomial
#'   outcomes must already be symmetric so averaging does not create invalid
#'   non-integer or non-binary observations. Influence covariates W must also be
#'   symmetric. Default is FALSE.
#'
#' @param bipartite Logical or NULL. Indicates whether the network is
#'   bipartite (senders and receivers are distinct node sets). If NULL
#'   (the default), bipartite status is inferred from Y: non-square arrays
#'   (n1 != n2) are treated as bipartite. Set to TRUE explicitly for
#'   square arrays where senders and receivers are nonetheless distinct
#'   populations. Setting FALSE on a non-square Y raises an error.
#'   Bipartite networks require \code{fix_receiver = TRUE} unless a separate
#'   \code{W_recv} is supplied (see below).
#'
#' @param W_recv Optional receiver-side influence covariate array
#'   (n2 x n2 x p2) for a \strong{full-bilinear bipartite} fit. When supplied,
#'   the sender influence \code{A} is built from \code{W} (n1 x n1 x p) and the
#'   receiver influence \code{B} from \code{W_recv}, fitting the complete
#'   \eqn{A X B'} model for two-mode data via alternating GLM (rather than
#'   collapsing to \code{B = I}). \code{alpha_1 = 1} pins the scale; all
#'   \code{beta} are free. Works for both rectangular (n1 != n2) and square
#'   two-mode networks. This path uses a pure-R alternating-GLM estimator (it
#'   does not call the C++ likelihood kernel), so it is slower than the square
#'   one-mode path for large networks; it draws random restarts (see
#'   \code{n_restarts} in \code{...}, default 5) and respects \code{seed}.
#'   Analytic standard errors are not available; use \code{\link{boot_sir}} with
#'   \code{type = "dyad"}. Default NULL.
#'
#' @param kron_mode Logical. \strong{Not yet implemented} (reserved for a future
#'   release): would estimate an unconstrained p x p coefficient matrix C instead
#'   of the rank-at-most-one alpha beta' factorization. Setting
#'   \code{kron_mode = TRUE} currently raises an error. Default is FALSE.
#'
#' @param seed Optional integer. If supplied, sets the random seed before
#'   fitting so that runs are reproducible (the estimators use random starting
#'   values). The global RNG state is restored on exit. Default NULL.
#'
#' @param ... Additional arguments passed to the fitting functions:
#'   \itemize{
#'     \item \code{trace}: Logical or integer controlling output verbosity.
#'     \item \code{tol}: Convergence tolerance for ALS (default 1e-8).
#'     \item \code{max_iter}: Maximum ALS iterations (default 100).
#'     \item \code{n_restarts}: Random restarts for the full-bilinear bipartite
#'       (\code{W_recv}) estimator (default 5); the lowest-deviance fit is kept.
#'   }
#'
#' @return An object of class \code{"sir"} with the following components:
#'   \describe{
#'     \item{summ}{Data frame of parameter estimates with columns \code{coef},
#'       \code{se} (classical SE), \code{rse} (robust/sandwich SE),
#'       \code{t_se} (z-statistic using classical SE),
#'       \code{t_rse} (z-statistic using robust SE). Row names identify each
#'       parameter.}
#'     \item{A}{Sender influence matrix. For static W: n1 x n1 matrix.
#'       For dynamic (4D) W: n1 x n1 x T array. Off-diagonal entry A[i,k]
#'       measures how much node k's behavior (via X) shapes node i's outgoing
#'       ties. Diagonal is set to zero for one-mode square fits; for a
#'       full-bilinear bipartite fit (\code{W_recv}) A is the dense
#'       \code{sum_k alpha_k W[,,k]} with all sender-side entries retained.}
#'     \item{B}{Receiver influence matrix. Identity when
#'       \code{fix_receiver = TRUE}; \code{n1 x n1} (same shape as A) for the
#'       square model; \code{n2 x n2} built from \code{W_recv} for a
#'       full-bilinear bipartite fit.}
#'     \item{tab}{Numeric vector of all estimated parameters in order:
#'       [theta_1, ..., theta_q, alpha_2, ..., alpha_p, beta_1, ..., beta_p].
#'       When \code{fix_receiver = TRUE}: [theta_1, ..., theta_q, alpha_1, ..., alpha_p].
#'       For a full-bilinear bipartite fit:
#'       [theta_1, ..., theta_q, alpha_2, ..., alpha_p, beta_1, ..., beta_p2].}
#'     \item{theta}{Coefficients for exogenous covariates Z (length q).}
#'     \item{alpha}{Full alpha vector including the fixed alpha_1 = 1 (length p).
#'       When \code{fix_receiver = TRUE}, all alpha are free.}
#'     \item{beta}{Coefficients for receiver influence covariates (length p, or
#'       length p2 for a full-bilinear bipartite fit from \code{W_recv}).
#'       Empty when \code{fix_receiver = TRUE}.}
#'     \item{p2}{Number of receiver-side influence covariates (only present for a
#'       full-bilinear bipartite fit).}
#'     \item{ll}{Log-likelihood at convergence.}
#'     \item{family}{The distribution family used (\code{"poisson"}, \code{"normal"},
#'       or \code{"binomial"}).}
#'     \item{method}{The estimation method used (\code{"ALS"} or \code{"optim"}).}
#'     \item{p}{Number of influence covariates in W.}
#'     \item{q}{Number of exogenous covariates in Z.}
#'     \item{m}{Number of sender nodes (same as n1).}
#'     \item{n1}{Number of sender (row) nodes.}
#'     \item{n2}{Number of receiver (column) nodes.}
#'     \item{bipartite}{Logical, TRUE if the network is treated as bipartite
#'       (non-square \code{Y}, \code{bipartite = TRUE}, or a \code{W_recv} fit).}
#'     \item{full_bilinear}{Logical, TRUE for a full-bilinear bipartite fit
#'       (\code{W_recv} supplied); such fits have no analytic SEs.}
#'     \item{n_periods}{Number of time periods.}
#'     \item{nobs}{Number of non-missing observations used in estimation.}
#'     \item{fitted.values}{Array (n1 x n2 x T) of fitted values on the response
#'       scale (counts for Poisson, probabilities for binomial, means for normal).}
#'     \item{residuals}{List with three components: \code{response} (Y - fitted),
#'       \code{pearson} (standardized by variance function), and
#'       \code{deviance} (signed square root of deviance contributions).}
#'     \item{vcov}{Variance-covariance matrix of parameters from the Hessian
#'       (classical SEs). NULL if \code{calc_se = FALSE}.}
#'     \item{vcov_robust}{Sandwich (robust) variance-covariance matrix.
#'       NULL if \code{calc_se = FALSE} or computation failed.}
#'     \item{Y}{The outcome array as used in fitting (with NAs from symmetric
#'       masking or Z missingness applied).}
#'     \item{W}{The influence covariate array.}
#'     \item{X}{The network state array (NAs replaced with 0).}
#'     \item{Z}{The exogenous covariate array (converted to 4D if 3D input).}
#'     \item{fix_receiver}{Logical, whether receiver effects were fixed.}
#'     \item{symmetric}{Logical, whether the network was treated as undirected.}
#'     \item{kron_mode}{Logical, whether Kronecker mode was used.}
#'     \item{iterations}{Number of iterations until convergence.}
#'     \item{history}{List with matrices ALPHA, BETA, THETA, DEV tracking
#'       parameter trajectories across iterations (useful for convergence
#'       diagnostics).}
#'     \item{convergence}{Logical, TRUE if the algorithm converged.}
#'     \item{call}{The matched function call.}
#'     \item{sigma2}{Estimated error variance (only for \code{family = "normal"}).}
#'   }
#'   
#' @references
#' Minhas, S. & Hoff, P. D. (2025). Social Influence Regression.
#' Political Analysis.
#' 
#' @examples
#' \donttest{
#' set.seed(123)
#' m <- 8; T_len <- 5; p <- 2
#' Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
#' X <- array(0, dim = c(m, m, T_len))
#' X[,,2:T_len] <- Y[,,1:(T_len - 1)]
#' W <- array(rnorm(m * m * p), dim = c(m, m, p))
#' model <- sir(Y = Y, W = W, X = X, family = "poisson",
#'              method = "ALS", calc_se = FALSE, max_iter = 10)
#' print(model)
#' coef(model)
#' }
#' 
#' @importFrom MASS ginv
#' @export
sir <- function(Y, W=NULL, X=NULL, Z=NULL, family, method="ALS", calc_se=TRUE,
				fix_receiver=FALSE, symmetric=FALSE, bipartite=NULL,
				W_recv=NULL, kron_mode=FALSE, seed=NULL, ...) {

	if (!family %in% c("poisson", "normal", "binomial")) {
	  cli::cli_abort("family must be one of {.val poisson}, {.val normal}, {.val binomial} (got {.val {family}}).")
	}

	# kron_mode is documented but not yet implemented; fail early and clearly.
	if (isTRUE(kron_mode)) {
	  cli::cli_abort(c(
		  "{.arg kron_mode} is not yet implemented.",
		  "i" = "Use the default rank-at-most-one model ({.code kron_mode = FALSE}); the unconstrained-C estimator is planned for a future release."
	  ))
	}

	# route full-bilinear bipartite fits
	if (!is.null(W_recv)) {
		if (isTRUE(symmetric)) {
			cli::cli_abort("{.arg symmetric} = TRUE is not compatible with full-bilinear bipartite fits.")
		}
		if (!is.null(bipartite) && !isTRUE(bipartite)) {
			cli::cli_abort("{.arg W_recv} fits are bipartite/two-mode; do not set {.arg bipartite} = FALSE.")
		}
		if (is.null(W) || is.null(X)) {
			cli::cli_abort("Full-bilinear bipartite fits require both {.arg W} and {.arg X}.")
		}
		dimsY <- dim(Y)
		if (length(dimsY) != 3) {
			cli::cli_abort("Y must be a 3D array with dimensions (n1 x n2 x T).")
		}
		n1 <- dimsY[1]; n2 <- dimsY[2]; T_len <- dimsY[3]
		if (any(is.infinite(Y) | is.nan(Y), na.rm = TRUE)) {
			n_inf <- sum(is.infinite(Y))
			n_nan <- sum(is.nan(Y))
			cli::cli_abort("Y contains non-finite values ({.val {n_inf}} Inf, {.val {n_nan}} NaN). Remove or impute these before fitting.")
		}
		Y_obs <- Y[!is.na(Y)]
		if (family == "binomial" && any(Y_obs != 0 & Y_obs != 1)) {
			n_bad <- sum(Y_obs != 0 & Y_obs != 1)
			cli::cli_abort("Binomial family requires Bernoulli outcomes coded 0/1; {.val {n_bad}} non-NA value{?s} are not 0 or 1.")
		}
		if (family == "poisson") {
			if (any(Y_obs < 0)) cli::cli_abort("Poisson family requires non-negative Y values.")
			if (any(abs(Y_obs - round(Y_obs)) > 1e-8)) {
				n_bad <- sum(abs(Y_obs - round(Y_obs)) > 1e-8)
				cli::cli_abort("Poisson family requires integer count outcomes; {.val {n_bad}} non-NA value{?s} are not integer counts.")
			}
		}
		if (!is.array(X) || length(dim(X)) != 3 || any(dim(X) != dimsY)) {
			cli::cli_abort("Dimension mismatch: X must have the same dimensions as Y ({.val {paste(dimsY, collapse = ' x ')}}).")
		}
		if (any(is.infinite(X) | is.nan(X), na.rm = TRUE)) {
			n_inf <- sum(is.infinite(X))
			n_nan <- sum(is.nan(X))
			cli::cli_abort("X contains non-finite values ({.val {n_inf}} Inf, {.val {n_nan}} NaN). Remove or impute these before fitting.")
		}
			if (!is.array(W) || length(dim(W)) != 3 || dim(W)[1] != n1 || dim(W)[2] != n1) {
				cli::cli_abort("{.arg W} must be {.val {n1}} x {.val {n1}} x p (sender covariates).")
			}
			if (!is.array(W_recv) || length(dim(W_recv)) != 3 || dim(W_recv)[1] != n2 || dim(W_recv)[2] != n2) {
				cli::cli_abort("{.arg W_recv} must be {.val {n2}} x {.val {n2}} x p2 (receiver covariates).")
			}
			if (dim(W)[3] < 1L) {
				cli::cli_abort("Full-bilinear bipartite fits require at least one sender-side influence covariate in {.arg W}.")
			}
			if (dim(W_recv)[3] < 1L) {
				cli::cli_abort("Full-bilinear bipartite fits require at least one receiver-side influence covariate in {.arg W_recv}.")
			}
		if (any(is.infinite(W) | is.nan(W), na.rm = TRUE)) {
			n_inf <- sum(is.infinite(W))
			n_nan <- sum(is.nan(W))
			cli::cli_abort("W contains non-finite values ({.val {n_inf}} Inf, {.val {n_nan}} NaN). Remove or impute these before fitting.")
		}
		if (any(is.infinite(W_recv) | is.nan(W_recv), na.rm = TRUE)) {
			n_inf <- sum(is.infinite(W_recv))
			n_nan <- sum(is.nan(W_recv))
			cli::cli_abort("W_recv contains non-finite values ({.val {n_inf}} Inf, {.val {n_nan}} NaN). Remove or impute these before fitting.")
		}
		if (!is.null(Z)) {
			dimsZ <- dim(Z)
			if (length(dimsZ) == 3) {
				if (!all(dimsZ == dimsY)) {
					cli::cli_abort("3D Z array must have dimensions matching Y ({.val {paste(dimsY, collapse = ' x ')}}), but got {.val {paste(dimsZ, collapse = ' x ')}}.")
				}
				Z <- array(Z, dim = c(n1, n2, 1, T_len),
						   dimnames = list(dimnames(Z)[[1]], dimnames(Z)[[2]], "Z1", dimnames(Z)[[3]]))
			} else if (length(dimsZ) == 4) {
				if (!all(dimsZ[c(1, 2, 4)] == dimsY)) {
					cli::cli_abort("4D Z array must have dimensions ({.val {n1}} x {.val {n2}} x q x {.val {T_len}}), but got {.val {paste(dimsZ, collapse = ' x ')}}.")
				}
			} else {
				cli::cli_abort("Z must be a 3D array (n1 x n2 x T) or 4D array (n1 x n2 x q x T), not {.val {length(dimsZ)}}D.")
			}
			if (any(is.infinite(Z) | is.nan(Z), na.rm = TRUE)) {
				n_inf <- sum(is.infinite(Z))
				n_nan <- sum(is.nan(Z))
				cli::cli_abort("Z contains non-finite values ({.val {n_inf}} Inf, {.val {n_nan}} NaN). Remove or impute these before fitting.")
			}
			if (anyNA(Z)) {
				z_na_mask <- apply(Z, c(1, 2, 4), function(x) any(is.na(x)))
				n_dyads_dropped <- sum(z_na_mask & !is.na(Y))
				if (n_dyads_dropped > 0) {
					cli::cli_inform("Z contains NAs: excluding {.val {n_dyads_dropped}} dyad-observations from likelihood (setting Y to NA where Z is missing).")
					Y[z_na_mask] <- NA
				}
				Z[is.na(Z)] <- 0
			}
		}
		if (anyNA(X)) X[is.na(X)] <- 0
		if (anyNA(W)) W[is.na(W)] <- 0
		if (anyNA(W_recv)) W_recv[is.na(W_recv)] <- 0

		y_send <- rownames(Y)
		y_recv <- colnames(Y)
		if (!is.null(y_send) || !is.null(y_recv)) {
			align_axis <- function(arr, axis, ref, ref_label, arr_label) {
				cur <- dimnames(arr)[[axis]]
				if (is.null(cur) || is.null(ref)) return(list(arr = arr, moved = FALSE))
				if (!setequal(cur, ref)) {
					cli::cli_abort(c(
						"Node names on {.arg {arr_label}} (dimension {axis}) do not match {.arg {ref_label}}.",
						"i" = "The name SETS differ, so the arrays describe different nodes; align them before fitting."
					))
				}
				if (identical(cur, ref)) return(list(arr = arr, moved = FALSE))
				idx <- match(ref, cur)
				d <- length(dim(arr))
				args <- rep(list(quote(expr = )), d)
				args[[axis]] <- idx
				list(arr = do.call(`[`, c(list(arr), args, list(drop = FALSE))), moved = TRUE)
			}
			reordered <- character(0)
			rX <- align_axis(X, 1, y_send, "Y", "X"); X <- rX$arr
			cX <- align_axis(X, 2, y_recv, "Y", "X"); X <- cX$arr
			if (rX$moved || cX$moved) reordered <- c(reordered, "X")
			if (!is.null(Z)) {
				rZ <- align_axis(Z, 1, y_send, "Y", "Z"); Z <- rZ$arr
				cZ <- align_axis(Z, 2, y_recv, "Y", "Z"); Z <- cZ$arr
				if (rZ$moved || cZ$moved) reordered <- c(reordered, "Z")
			}
			rW <- align_axis(W, 1, y_send, "Y", "W"); W <- rW$arr
			cW <- align_axis(W, 2, y_send, "Y", "W"); W <- cW$arr
			if (rW$moved || cW$moved) reordered <- c(reordered, "W")
			rWr <- align_axis(W_recv, 1, y_recv, "Y", "W_recv"); W_recv <- rWr$arr
			cWr <- align_axis(W_recv, 2, y_recv, "Y", "W_recv"); W_recv <- cWr$arr
			if (rWr$moved || cWr$moved) reordered <- c(reordered, "W_recv")
			if (length(reordered) > 0) {
				cli::cli_inform("Reordered {.arg {reordered}} node dimension{?s} to match {.arg Y} (dimnames were present but in a different order).")
			}
		}

		if (isTRUE(calc_se)) {
			cli::cli_inform(c(
				"Analytic standard errors are unavailable for full-bilinear bipartite fits.",
				"i" = "Use {.fn boot_sir} with {.code type = \"dyad\"} for inference."
			))
		}
		# set a local random seed
		if (!is.null(seed)) {
			if (!is.numeric(seed) || length(seed) != 1) {
				cli::cli_abort("{.arg seed} must be a single number or NULL.")
			}
			if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
				old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
				on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
			} else {
				on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
					rm(".Random.seed", envir = globalenv()), add = TRUE)
			}
			set.seed(seed)
		}
		dots <- list(...)
		res <- .sir_bipartite(Y, W, W_recv, X, Z, family,
			max_iter = if (!is.null(dots$max_iter)) dots$max_iter else 50L,
			tol = if (!is.null(dots$tol)) dots$tol else 1e-6,
			n_restarts = if (!is.null(dots$n_restarts)) dots$n_restarts else 5L,
			trace = isTRUE(dots$trace))
		res$seed <- seed
		res$call <- match.call()
		return(res)
	}

	# normalize the estimation method
	method <- switch(tolower(method[1]),
	als = "ALS",
	optim = ,
	bfgs = "optim",
	cli::cli_abort(c(
	  "{.arg method} must be one of {.val ALS} or {.val optim}.",
	  "i" = "{.val bfgs} is accepted as an alias for {.val optim}; matching is case-insensitive."
	))
	)

	# set a local random seed
	if (!is.null(seed)) {
	  if (!is.numeric(seed) || length(seed) != 1) {
		  cli::cli_abort("{.arg seed} must be a single number or NULL.")
	  }
	  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
		  old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
		  on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
	  } else {
		  on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
			  rm(".Random.seed", envir = globalenv()), add = TRUE)
	  }
	  set.seed(seed)
	}

	# input validation and dimension checks
	dimsY <- dim(Y)
	if (length(dimsY) != 3) {
	  cli::cli_abort("Y must be a 3D array with dimensions (n1 x n2 x T).")
	}

	# check for non-finite values
	if (any(is.infinite(Y) | is.nan(Y), na.rm = TRUE)) {
	  n_inf <- sum(is.infinite(Y))
	  n_nan <- sum(is.nan(Y))
	  cli::cli_abort("Y contains non-finite values ({.val {n_inf}} Inf, {.val {n_nan}} NaN). Remove or impute these before fitting.")
	}

	# validate Y for binomial family
	if (family == "binomial") {
	  Y_obs <- Y[!is.na(Y)]
	  if (any(Y_obs != 0 & Y_obs != 1)) {
		  n_bad <- sum(Y_obs != 0 & Y_obs != 1)
		  cli::cli_abort("Binomial family requires Bernoulli outcomes coded 0/1; {.val {n_bad}} non-NA value{?s} are not 0 or 1.")
	  }
	}

	# validate Y for poisson family
	if (family == "poisson") {
	  Y_obs <- Y[!is.na(Y)]
	  if (any(Y_obs < 0)) {
		  cli::cli_abort("Poisson family requires non-negative Y values.")
	  }
	  if (any(abs(Y_obs - round(Y_obs)) > 1e-8)) {
		  n_bad <- sum(abs(Y_obs - round(Y_obs)) > 1e-8)
		  cli::cli_abort("Poisson family requires integer count outcomes; {.val {n_bad}} non-NA value{?s} are not integer counts.")
	  }
	}

	n1 <- dimsY[1]
	n2 <- dimsY[2]
	T_len <- dimsY[3]

	# bipartite detection: user-supplied value takes precedence, otherwise infer from dimensions
		if (is.null(bipartite)) {
		  bipartite <- (n1 != n2)
		  if (!bipartite &&
			  !is.null(dimnames(Y)[[1]]) && !is.null(dimnames(Y)[[2]]) &&
			  !identical(dimnames(Y)[[1]], dimnames(Y)[[2]])) {
			  cli::cli_abort(c(
				  "Square {.arg Y} has different row and column labels but {.arg bipartite} was not set.",
				  "i" = "Use {.code bipartite = TRUE} when rows and columns are distinct actor sets, even if their counts are equal."
			  ))
		  }
		} else {
	  if (!is.logical(bipartite) || length(bipartite) != 1) {
		  cli::cli_abort("{.arg bipartite} must be TRUE, FALSE, or NULL (auto-detect).")
	  }
	  if (bipartite && n1 == n2) {
		  cli::cli_inform("Treating square {.val {n1}} x {.val {n2}} network as bipartite (user-specified).")
	  }
	  if (!bipartite && n1 != n2) {
		  cli::cli_abort("Y is non-square ({.val {n1}} x {.val {n2}}) but {.arg bipartite} = FALSE. Non-square arrays are necessarily bipartite.")
	  }
	}

	if (bipartite && symmetric) {
	  cli::cli_abort("{.arg symmetric} = TRUE is not compatible with bipartite (non-square) networks.")
	}

	one_mode_square <- !bipartite && n1 == n2
	if (one_mode_square) {
	  Y <- set_square_diagonal(Y, NA)
	}

	if (bipartite) {
	  if (!fix_receiver) {
		  cli::cli_inform(c(
			  "Bipartite network with sender-side {.arg W} only: fixing {.arg fix_receiver} = TRUE (B = I).",
			  "i" = "To estimate the full bilinear A*X*B' with a separate receiver-side influence, supply {.arg W_recv} (an n2 x n2 x p2 array)."
		  ))
		  fix_receiver <- TRUE
	  }
	  if (method == "optim") {
		  cli::cli_inform("Bipartite networks use ALS method. Switching from {.val optim} to {.val ALS}.")
		  method <- "ALS"
	  }
	}

	m <- n1  # sender dimension (= n2 for square networks)

	# warn when T=1 and influence structure is requested
	if (T_len == 1 && !is.null(W) && !all(dim(W)[3] == 0)) {
	  cli::cli_warn("Only one time period (T = 1). If X is all zeros the influence term A*X*B' contributes nothing. Consider whether the model is identified.")
	}

	# symmetric network handling
	if (symmetric) {
	  # check symmetry before masking the triangle
	  Y_sym_diff <- max(abs(Y - aperm(Y, c(2,1,3))), na.rm = TRUE)
	  if (Y_sym_diff > 1e-10) {
		  if (family %in% c("poisson", "binomial")) {
			  cli::cli_abort(c(
				  "For {.arg symmetric} = TRUE with {.val {family}} outcomes, {.arg Y} must already be symmetric.",
				  "i" = "Averaging asymmetric discrete outcomes can create invalid non-integer or non-binary values."
			  ))
		  }
		  cli::cli_inform("Symmetrizing Y: max|Y - Y'| = {.val {sprintf('%.4f', Y_sym_diff)}}. Averaging upper and lower triangles.")
		  Y <- (Y + aperm(Y, c(2,1,3))) / 2
	  }
	  # use only upper-triangle observations
	  for (t in 1:T_len) {
		  Y[,,t][lower.tri(Y[,,t])] <- NA
		  diag(Y[,,t]) <- NA
	  }
	  # symmetrize X if provided
	  if (!is.null(X)) {
		  X <- (X + aperm(X, c(2,1,3))) / 2
		  X[is.na(X)] <- 0
	  }
	  # undirected: only sender-side influence estimated
	  fix_receiver <- TRUE
	}

	# detect dynamic (4D) vs static (3D) W
	dynamic_W <- FALSE
	if (!is.null(W)) {
	  dimsW <- dim(W)
	  if (length(dimsW) == 4) {
		  # W is (m x m x p x T): time-varying influence covariates
		  dynamic_W <- TRUE
		  p <- dimsW[3]
		  if (dimsW[4] != T_len) {
			  cli::cli_abort("4D W has {.val {dimsW[4]}} time slices but Y has {.val {T_len}}. They must match.")
		  }
		  if (bipartite) {
			  if (dimsW[1] != n1 || dimsW[2] != n1) {
				  cli::cli_abort("For bipartite networks, 4D W must be {.val {n1}} x {.val {n1}} x p x T (sender covariates), but got {.val {dimsW[1]}} x {.val {dimsW[2]}} x {.val {dimsW[3]}} x {.val {dimsW[4]}}.")
			  }
		  } else if (any(dimsW[1:2] != dimsY[1:2])) {
			  cli::cli_abort("Dimension mismatch: Y is {.val {dimsY[1]}} x {.val {dimsY[2]}} but W is {.val {dimsW[1]}} x {.val {dimsW[2]}}.")
		  }
	  } else if (length(dimsW) == 3) {
		  p <- dimsW[3]
	  } else {
		  cli::cli_abort("W must be a 3D array (m x m x p) or 4D array (m x m x p x T), not {.val {length(dimsW)}}D.")
	  }

	  if (is.null(X)) {
		  cli::cli_abort("X must be provided when W is specified. X typically contains lagged Y values.")
	  }
	  if (bipartite) {
		  # bipartite: W must be n1 x n1 x p (sender covariates)
		  if (dim(W)[1] != n1 || dim(W)[2] != n1) {
			  cli::cli_abort("For bipartite networks, W must be {.val {n1}} x {.val {n1}} x p (sender covariates), but got {.val {dim(W)[1]}} x {.val {dim(W)[2]}} x {.val {dim(W)[3]}}.")
		  }
	  } else if (!dynamic_W) {
		  if (any(dimsY[1:2] != dim(W)[1:2])) {
			  cli::cli_abort("Dimension mismatch: Y is {.val {dimsY[1]}} x {.val {dimsY[2]}} but W is {.val {dim(W)[1]}} x {.val {dim(W)[2]}}. The first two dimensions must match.")
		  }
	  }
	  if (any(is.infinite(X) | is.nan(X), na.rm = TRUE)) {
		  n_inf <- sum(is.infinite(X))
		  n_nan <- sum(is.nan(X))
		  cli::cli_abort("X contains non-finite values ({.val {n_inf}} Inf, {.val {n_nan}} NaN). Remove or impute these before fitting.")
	  }
	  if (any(dimsY != dim(X))) {
		  cli::cli_abort("Dimension mismatch: Y is {.val {paste(dimsY, collapse = ' x ')}} but X is {.val {paste(dim(X), collapse = ' x ')}}. They must have the same dimensions.")
	  }
	  # replace missing influence states with zero
	  n_na_X <- sum(is.na(X))
	  if (n_na_X > 0) {
		  n_na_X_offdiag <- n_na_X
		  if (one_mode_square) {
			  for (t in seq_len(T_len)) n_na_X_offdiag <- n_na_X_offdiag - sum(is.na(diag(X[,,t])))
		  }
		  if (n_na_X_offdiag > 0) {
			  cli::cli_inform("Replacing {.val {n_na_X_offdiag}} NA values in X with 0 (missing observations carry zero influence).")
		  }
		  X[is.na(X)] <- 0
	  }
		  if (one_mode_square) {
			  X <- set_square_diagonal(X, 0)
		  }

	  # replace missing influence covariates with zero
	  n_na_W <- sum(is.na(W))
	  if (n_na_W > 0) {
		  if (!dynamic_W && one_mode_square) {
			  n_na_W_offdiag <- n_na_W
			  for (k in seq_len(p)) n_na_W_offdiag <- n_na_W_offdiag - sum(is.na(diag(W[,,k])))
		  } else {
			  n_na_W_offdiag <- n_na_W
		  }
		  if (n_na_W_offdiag > 0) {
			  cli::cli_inform("Replacing {.val {n_na_W_offdiag}} NA values in W with 0 (missing covariates excluded from influence structure).")
		  }
		  W[is.na(W)] <- 0
	  }
		  if (one_mode_square) {
			  W <- set_square_diagonal(W, 0)
		  }
		  if (symmetric) {
			  max_w_diff <- 0
			  if (dynamic_W) {
				  for (tt in seq_len(T_len)) {
					  for (rr in seq_len(p)) {
						  Wrt <- W[, , rr, tt]
						  max_w_diff <- max(max_w_diff, max(abs(Wrt - t(Wrt)), na.rm = TRUE))
					  }
				  }
			  } else {
				  for (rr in seq_len(p)) {
					  Wr <- W[, , rr]
					  max_w_diff <- max(max_w_diff, max(abs(Wr - t(Wr)), na.rm = TRUE))
				  }
			  }
			  if (is.finite(max_w_diff) && max_w_diff > 1e-10) {
				  cli::cli_abort(c(
					  "For {.arg symmetric} = TRUE, every influence-covariate slice in {.arg W} must be symmetric.",
					  "i" = "Symmetrize W before fitting, for example {.code W[, , r] <- (W[, , r] + t(W[, , r])) / 2}."
				  ))
			  }
		  }

	  # warn about W collinearity (can cause singular hessian)
	  # for dynamic W, check the first time slice as representative.
	  # gate on p >= 2 (the minimal common model) so collinear 2-covariate W is
	  # also flagged -- the SVD condition number needs at least 2 columns.
	  if (p >= 2 && !bipartite) local({
		  W_check <- if (dynamic_W) W[,,,1] else W
		  W_vecs <- matrix(NA, m * (m - 1), p)
		  for (k in 1:p) {
			  wk <- W_check[,,k]
			  W_vecs[, k] <- wk[row(wk) != col(wk)]
		  }
		  # only scale columns with nonzero variance
		  W_sd <- apply(W_vecs, 2, sd, na.rm = TRUE)
		  W_vecs_use <- W_vecs[, W_sd > 0, drop = FALSE]
		  if (ncol(W_vecs_use) > 1) {
			  W_svd <- svd(scale(W_vecs_use))$d
			  W_cond <- W_svd[1] / max(W_svd[length(W_svd)], .Machine$double.eps)
			  if (W_cond > 100) {
				  cli::cli_warn("W covariates have high collinearity (condition number: {.val {sprintf('%.1f', W_cond)}}). This may cause singular Hessian and unreliable SEs.")
			  }
		  }
	  })
	} else {
	  p <- 0
	  # dummy empty arrays for C++ when p=0
	  W <- array(0, dim=c(n1, n1, 0))
	  # X still supplies dimensions for eta_tab even with p=0.
	  if (is.null(X)) {
		  X <- array(0, dim=c(n1, n2, T_len))
	  } else if (!is.array(X) || length(dim(X)) != 3L || any(dim(X) != dimsY)) {
		  cli::cli_abort("X was supplied without W, but its dimensions do not match Y ({.val {paste(dimsY, collapse = ' x ')}}).")
	  }
	}

	# handle Z dimensions (allow 3D or 4D)
	if (!is.null(Z)) {
	  dimsZ <- dim(Z)
	  if (length(dimsZ) == 3) {
		  if (all(dimsZ == dimsY)) {
			  # store one covariate with an explicit covariate axis
			  Z <- array(Z, dim=c(n1, n2, 1, T_len))
			  q <- 1
		  } else {
			  cli::cli_abort("3D Z array must have dimensions matching Y ({.val {paste(dimsY, collapse = ' x ')}}), but got {.val {paste(dimsZ, collapse = ' x ')}}.")
		  }
	  } else if (length(dimsZ) == 4) {
		  if (all(dimsZ[c(1,2,4)] == dimsY)) {
			  q <- dimsZ[3]
		  } else {
			  cli::cli_abort("4D Z array must have dimensions ({.val {n1}} x {.val {n2}} x q x {.val {T_len}}), but got {.val {paste(dimsZ, collapse = ' x ')}}.")
		  }
	  } else {
		  cli::cli_abort("Z must be a 3D array (m x m x T) or 4D array (m x m x q x T), not {.val {length(dimsZ)}}D.")
	  }
	  if (symmetric) {
		  max_z_diff <- 0
		  for (tt in seq_len(T_len)) {
			  for (rr in seq_len(q)) {
				  Zrt <- Z[, , rr, tt]
				  if (!all(is.na(Zrt) == t(is.na(Zrt)))) {
					  cli::cli_abort(c(
						  "For {.arg symmetric} = TRUE, every direct-covariate slice in {.arg Z} must be symmetric.",
						  "i" = "Symmetrize Z before fitting so each unordered pair has one value."
					  ))
				  }
				  finite_pair <- !is.na(Zrt) & !is.na(t(Zrt))
				  slice_diff <- if (any(finite_pair)) {
					  max(abs(Zrt[finite_pair] - t(Zrt)[finite_pair]), na.rm = TRUE)
				  } else {
					  0
				  }
				  max_z_diff <- max(max_z_diff, slice_diff)
			  }
		  }
		  if (is.finite(max_z_diff) && max_z_diff > 1e-10) {
			  cli::cli_abort(c(
				  "For {.arg symmetric} = TRUE, every direct-covariate slice in {.arg Z} must be symmetric.",
				  "i" = "Symmetrize Z before fitting, for example {.code Z[, , r, t] <- (Z[, , r, t] + t(Z[, , r, t])) / 2}."
			  ))
		  }
	  }
	} else {
	  q <- 0
	  # Z remains NULL
	}

	# handle NAs in Z: propagate to Y (excludes dyad from likelihood), then zero out
	if (!is.null(Z) && anyNA(Z)) {
	  # dyad is missing if any of its q covariates is NA
	  z_na_mask <- apply(Z, c(1, 2, 4), function(x) any(is.na(x)))
	  n_dyads_dropped <- sum(z_na_mask & !is.na(Y))
	  if (n_dyads_dropped > 0) {
		  cli::cli_inform("Z contains NAs: excluding {.val {n_dyads_dropped}} dyad-observations from likelihood (setting Y to NA where Z is missing).")
		  Y[z_na_mask] <- NA
	  }
	  Z[is.na(Z)] <- 0
	}

	# align Y/W/X/Z by node dimnames when present; reorder W/X/Z to Y when the
	# name sets match but order differs, and abort when the sets differ
	y_send <- rownames(Y)   # sender (row) names; NULL if absent
	y_recv <- colnames(Y)   # receiver (col) names; NULL if absent
	if (!is.null(y_send) || !is.null(y_recv)) {
	  # reorder one array axis to match a reference name vector; returns the
	  # (possibly reordered) array and a flag noting whether anything moved.
	  align_axis <- function(arr, axis, ref, ref_label, arr_label) {
		  cur <- dimnames(arr)[[axis]]
		  if (is.null(cur) || is.null(ref)) return(list(arr = arr, moved = FALSE))
		  if (!setequal(cur, ref)) {
			  cli::cli_abort(c(
				  "Node names on {.arg {arr_label}} (dimension {axis}) do not match {.arg {ref_label}}.",
				  "i" = "The name SETS differ, so the arrays describe different nodes; align them before fitting.",
				  "i" = "In {.arg {arr_label}} only: {.val {setdiff(cur, ref)}}.",
				  "i" = "In {.arg {ref_label}} only: {.val {setdiff(ref, cur)}}."
			  ))
		  }
		  if (identical(cur, ref)) return(list(arr = arr, moved = FALSE))
		  # same set, different order: reorder this axis to match ref
		  idx <- match(ref, cur)
		  d <- length(dim(arr))
		  args <- rep(list(quote(expr = )), d)
		  args[[axis]] <- idx
		  arr2 <- do.call(`[`, c(list(arr), args, list(drop = FALSE)))
		  list(arr = arr2, moved = TRUE)
	  }

	  reordered <- character(0)

	  # X: dims 1,2 = senders, receivers (same axes as Y)
	  if (!is.null(X)) {
		  rX <- align_axis(X, 1, y_send, "Y", "X"); X <- rX$arr
		  cX <- align_axis(X, 2, y_recv, "Y", "X"); X <- cX$arr
		  if (rX$moved || cX$moved) reordered <- c(reordered, "X")
	  }

	  # Z (4D here): dims 1,2 = senders, receivers
	  if (!is.null(Z)) {
		  rZ <- align_axis(Z, 1, y_send, "Y", "Z"); Z <- rZ$arr
		  cZ <- align_axis(Z, 2, y_recv, "Y", "Z"); Z <- cZ$arr
		  if (rZ$moved || cZ$moved) reordered <- c(reordered, "Z")
	  }

	  # W: bipartite/static-square W is sender-by-sender (both axes = Y's row
	  # order); for square directed W dims 1,2 are sender,receiver. dynamic 4D W
	  # carries nodes on dims 1,2. the p=0 dummy W has no dimnames so is untouched.
	  if (p > 0) {
		  if (bipartite) {
			  rW <- align_axis(W, 1, y_send, "Y", "W"); W <- rW$arr
			  cW <- align_axis(W, 2, y_send, "Y", "W"); W <- cW$arr
			  if (rW$moved || cW$moved) reordered <- c(reordered, "W")
		  } else {
			  rW <- align_axis(W, 1, y_send, "Y", "W"); W <- rW$arr
			  cW <- align_axis(W, 2, y_recv, "Y", "W"); W <- cW$arr
			  if (rW$moved || cW$moved) reordered <- c(reordered, "W")
		  }
	  }

	  if (length(reordered) > 0) {
		  cli::cli_inform("Reordered {.arg {reordered}} node dimension{?s} to match {.arg Y} (dimnames were present but in a different order).")
	  }
	}

	# report overall NA pattern in Y (off-diagonal only, diagonals always excluded)
	if (one_mode_square) {
	  n_na_offdiag <- 0L
	  for (t in seq_len(T_len)) {
		  Yt <- Y[,,t]
		  diag(Yt) <- 0  # don't count diagonal
		  n_na_offdiag <- n_na_offdiag + sum(is.na(Yt))
	  }
	  if (n_na_offdiag > 0) {
		  n_offdiag <- n1 * (n1 - 1) * T_len
		  pct_miss <- sprintf("%.1f", 100 * n_na_offdiag / n_offdiag)
		  cli::cli_inform("Fitting on {.val {n_offdiag - n_na_offdiag}} of {.val {n_offdiag}} off-diagonal dyad-observations ({pct_miss}% missing).")
	  }
	} else {
	  n_na_Y <- sum(is.na(Y))
	  if (n_na_Y > 0) {
		  n_total <- prod(dim(Y))
		  pct_miss <- sprintf("%.1f", 100 * n_na_Y / n_total)
		  cli::cli_inform("Fitting on {.val {n_total - n_na_Y}} of {.val {n_total}} dyad-observations ({pct_miss}% missing).")
	  }
	}

	# fitting
	if (fix_receiver && method == "optim") {
	cli::cli_warn("{.arg fix_receiver} = TRUE forces ALS method (single GLM step).")
	method <- "ALS"
	}

	if (dynamic_W && method == "optim") {
	  cli::cli_inform("Dynamic (4D) W uses ALS method. Switching from {.val optim} to {.val ALS}.")
	  method <- "ALS"
	}

	if (method == "ALS") {
	mod <- sir_alsfit(Y, W, X, Z, family,
					  fix_receiver=fix_receiver, kron_mode=kron_mode,
					  dynamic_W=dynamic_W, ...)
	} else if (method == "optim") {
	mod <- sir_optfit(Y, W, X, Z, family, ...)
	} else {
	cli::cli_abort("method must be {.val ALS} or {.val optim} (got {.val {method}}).")
	}

	tab <- mod$tab
	fit_converged <- if (method == "optim") {
		isTRUE(mod$convergence == 0)
	} else {
		isTRUE(mod$converged)
	}

	# reconstruct final alpha and beta
	if (fix_receiver && p > 0) {
	  # parse the fixed-receiver parameter vector
	  if (q > 0) theta <- tab[1:q] else theta <- numeric(0)
	  alpha <- tab[(q+1):(q+p)]
	  beta <- numeric(0)
	} else {
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
			  alpha <- c(1, a)
		  } else {
			  a <- numeric(0)
			  alpha <- c(1)
		  }
		  b_start <- q + p
		  beta <- tab[b_start:length(tab)]
	  } else {
		  alpha <- numeric(0)
		  beta <- numeric(0)
	  }
	}


	# influence matrices A, B
	# sign ambiguity is resolved by alpha_1 = 1 constraint; no additional
	# sign flipping is needed here
		if (dynamic_W && p > 0) {
		  # dynamic W: A_t and B_t are time-varying (3D arrays)
		  A <- array(0, dim = c(n1, n1, T_len))
		  B <- array(0, dim = c(n2, n2, T_len))
		  for (t in seq_len(T_len)) {
			  W_t <- W[,,,t]
			  A[,,t] <- cpp_amprod_W_v(W_t, alpha)
		  if (fix_receiver) {
			  B[,,t] <- diag(n2)
			  } else {
				  B[,,t] <- cpp_amprod_W_v(W_t, beta)
			  }
			  if (!bipartite) diag(A[,,t]) <- 0
				  if (!fix_receiver) diag(B[,,t]) <- 0
			  }
		  send_nm <- dimnames(Y)[[1]]
		  recv_nm <- dimnames(Y)[[2]]
		  time_nm <- dimnames(Y)[[3]]
		  if (!is.null(send_nm)) dimnames(A) <- list(send_nm, send_nm, time_nm)
		  if (!is.null(recv_nm)) dimnames(B) <- list(recv_nm, recv_nm, time_nm)
		} else if (fix_receiver && p > 0) {
	  A <- cpp_amprod_W_v(W, alpha)
	  B <- diag(n2)
	  if (!bipartite) diag(A) <- 0

	  if (!is.null(rownames(Y))) {
		rownames(A) <- colnames(A) <- rownames(Y)
	  }
	  if (!is.null(colnames(Y))) {
		rownames(B) <- colnames(B) <- colnames(Y)
	  }
	} else if (p > 0) {
	  A <- cpp_amprod_W_v(W, alpha)
	  B <- cpp_amprod_W_v(W, beta)
	  diag(A) <- 0
	  diag(B) <- 0

	  if (!is.null(rownames(Y))) {
		rownames(A) <- colnames(A) <- rownames(Y)
		rownames(B) <- colnames(B) <- rownames(Y)
	  }
	} else {
	  A <- matrix(0, n1, n1)
	  B <- matrix(0, n2, n2)
	}


	# log-likelihood and variance estimation
		NLL <- mll_sir(tab, Y, W, X, Z, family,
					   fix_receiver = fix_receiver, bipartite = bipartite)
	ll <- -NLL

	# fitted values on response scale
	eta <- eta_tab(tab, W, X, Z, fix_receiver=fix_receiver)
	if (family == "poisson") {
	  fitted_values <- exp(eta)
	} else if (family == "binomial") {
	  prob <- 1 / (1 + exp(-eta))
	  prob[prob < 1e-15] <- 1e-15
	  prob[prob > 1 - 1e-15] <- 1 - 1e-15
	  fitted_values <- prob
	} else {
	  fitted_values <- eta
	}

		if (symmetric) {
		  for (t in seq_len(T_len)) {
			  upper <- fitted_values[, , t]
			  upper[lower.tri(upper, diag = TRUE)] <- NA
			  upper0 <- upper
			  upper0[is.na(upper0)] <- 0
			  mirrored <- upper0 + t(upper0)
			  diag(mirrored) <- NA
			  fitted_values[, , t] <- mirrored
		  }
		} else if (one_mode_square) {
		  fitted_values <- set_square_diagonal(fitted_values, NA)
		}

	sigma2_est <- 1.0
	sigma2_mle <- NULL

	if (family == "normal") {
	  # estimate gaussian variance from likelihood observations
	  N_obs_ll <- 0L
		  if (one_mode_square) {
			  for (t in seq_len(T_len)) {
				  Yt <- Y[,,t]
				  diag(Yt) <- NA
			  N_obs_ll <- N_obs_ll + sum(!is.na(Yt))
		  }
	  } else {
		  N_obs_ll <- sum(!is.na(Y))
	  }
	  df <- N_obs_ll - length(tab)
	  # guard against an unidentified normal fit (too few obs / collinear design)
	  if (df <= 0) {
		  cli::cli_abort(c(
			  "Non-positive residual degrees of freedom for the normal family (N = {N_obs_ll}, parameters = {length(tab)}).",
			  "i" = "The model is not identified with this much data; reduce p/q or supply more observations."
		  ))
	  }
	  RSS <- sum((Y - eta)^2, na.rm = TRUE)
	  sigma2_est <- RSS / df
	  if (!is.finite(sigma2_est) || sigma2_est <= 0) {
		  cli::cli_abort(c(
			  "Could not estimate a valid residual variance (sigma^2) for the normal family.",
			  "i" = "Estimated sigma^2 = {.val {sigma2_est}}; check for a perfect fit, too few observations, or collinear covariates."
		  ))
	  }
	  # Full Gaussian likelihood for IC/logLik uses the MLE scale RSS / N.
	  # Keep sigma2_est as the residual-variance estimate used for SE scaling.
	  sigma2_mle <- RSS / N_obs_ll
	  ll <- -0.5 * N_obs_ll * log(2 * pi * sigma2_mle) - 0.5 * RSS / sigma2_mle
	}


	# standard errors
	summ <- data.frame(coef=tab)
	vcov_mat <- NULL
	vcov_robust <- NULL
	# hessian-conditioning diagnostics, separate from solver convergence; default
	# NA so the fields always exist, populated in the se section below
	hessian_cond_num <- NA_real_
	hessian_min_eig <- NA_real_
	se_reliable <- NA
	if (calc_se && !fit_converged) {
		cli::cli_warn(c(
			"Model did not converge; standard errors are suppressed.",
			"i" = "Refit until convergence before reporting Wald inference."
		))
		se_reliable <- FALSE
		calc_se <- FALSE
	}

	if(calc_se){

	if (fix_receiver && p > 0 && !is.null(mod$glm_fit)) {
	  # use the fitted glm covariance
	  glm_fit <- mod$glm_fit

	  vcov_mat <- tryCatch(vcov(glm_fit), error = function(e) NULL)

	  # mark the glm covariance as reliable when finite
	  if (!is.null(vcov_mat)) {
		se <- sqrt(pmax(diag(vcov_mat), 0))
		se_reliable <- all(is.finite(diag(vcov_mat)))

		# hc0 sandwich rebuilt from the design and response the glm step used; the
		# score factor is the response residual (y - mu), bread is the glm vcov
		vcov_robust <- tryCatch({
		  Xd <- mod$glm_X
		  yv <- mod$glm_y
		  eta_fr <- as.numeric(Xd %*% tab)
		  mu_fr <- switch(family,
						  poisson  = exp(eta_fr),
						  binomial = 1 / (1 + exp(-eta_fr)),
						  normal   = eta_fr)
		  r_fr <- yv - mu_fr
		  # rows the GLM used = complete response (NA y dropped, incl. diagonal)
		  complete <- !is.na(r_fr) & is.finite(r_fr)
		  X_c <- Xd[complete, , drop = FALSE]
		  r_c <- r_fr[complete]
		  meat <- crossprod(X_c * r_c)
		  bread <- vcov_mat / ifelse(family == "normal", sigma2_est, 1)
		  bread %*% meat %*% bread
		}, error = function(e) NULL)

		rse <- tryCatch({
		  rse_diag <- diag(vcov_robust)
		  rse_diag[rse_diag < 0] <- NA
		  sqrt(rse_diag)
		}, error = function(e) rep(NA, length(se)))

		summ$se <- se
		summ$rse <- rse
		summ$t_se <- summ$coef / summ$se
		summ$t_rse <- summ$coef / summ$rse
	  } else {
		se_reliable <- FALSE
		summ$se <- NA
		summ$rse <- NA
		summ$t_se <- NA
		summ$t_rse <- NA
	  }

	} else if (dynamic_W) {
	  # dynamic W: use C++ backend with list-of-cubes representation
	  W_field <- prepare_W_field(W)
	  Z_list <- prepare_Z_list(Z)
	  gH <- cpp_mll_gH_dyn(tab, Y, W_field, X, Z_list, family)

	  H <- gH$hess
	  if (family == "normal") {
		  H <- H / sigma2_est
		  S <- gH$shess / (sigma2_est^2)
	  } else {
		  S <- gH$shess
	  }

	  H_eig <- eigen(H, symmetric = TRUE)
	  min_eig <- min(H_eig$values)
	  max_eig <- max(H_eig$values)
	  cond_num <- max_eig / max(max(min_eig, .Machine$double.eps), .Machine$double.eps)

	  # record hessian conditioning; se_reliable is false when ill-conditioned
	  hessian_min_eig <- min_eig
	  hessian_cond_num <- cond_num
	  se_reliable <- !(min_eig <= 0 || cond_num > 1e10)

	  if (min_eig <= 0 || cond_num > 1e10) {
		  cli::cli_warn("Hessian is ill-conditioned (min eigenvalue: {.val {sprintf('%.2e', min_eig)}}, condition number: {.val {sprintf('%.2e', cond_num)}}). Applying ridge regularization. Consider bootstrap SEs via {.fn boot_sir} for reliable inference.")
		  lambda <- max_eig * 1e-6
		  H_reg <- H + lambda * diag(nrow(H))
		  H_inv <- tryCatch(solve(H_reg), error = function(e) {
			  tryCatch(MASS::ginv(H_reg), error = function(e2) NULL)
		  })
	  } else {
		  H_inv <- tryCatch(solve(H), error = function(e) {
			  tryCatch(MASS::ginv(H), error = function(e2) NULL)
		  })
	  }

	  if (!is.null(H_inv)) {
		  se_diag <- diag(H_inv)
		  se_diag[se_diag < 0] <- NA
		  se <- sqrt(se_diag)
		  vcov_mat <- H_inv
		  vcov_robust <- tryCatch(H_inv %*% S %*% H_inv, error = function(e) NULL)
		  rse <- tryCatch({
			  rse_diag <- diag(vcov_robust)
			  rse_diag[rse_diag < 0] <- NA
			  sqrt(rse_diag)
		  }, error = function(e) rep(NA, length(se)))

		  summ$se <- se
		  summ$rse <- rse
		  summ$t_se <- summ$coef / summ$se
		  summ$t_rse <- summ$coef / summ$rse
	  } else {
		  se_reliable <- FALSE
		  summ$se <- NA; summ$rse <- NA; summ$t_se <- NA; summ$t_rse <- NA
	  }

	} else {
	  # standard path: use C++ backend for hessian
	  Z_list <- prepare_Z_list(Z)
	  gH <- cpp_mll_gH(tab, Y, W, X, Z_list, family)

	  H <- gH$hess

	  # adjust hessian/score for gaussian sigma^2
	  if (family == "normal") {
		  # H_true = H_cpp / sigma^2, S_true = S_cpp / sigma^4
		  H <- H / sigma2_est
		  # C++ assumes sigma=1, adjust here
		  S <- gH$shess / (sigma2_est^2)
	  } else {
		  S <- gH$shess
	  }


	  # check hessian invertibility via eigendecomposition
	  H_eig <- eigen(H, symmetric = TRUE)
	  min_eig <- min(H_eig$values)
	  max_eig <- max(H_eig$values)
	  cond_num <- max_eig / max(max(min_eig, .Machine$double.eps), .Machine$double.eps)

	  # record hessian conditioning; se_reliable is false when ill-conditioned
	  hessian_min_eig <- min_eig
	  hessian_cond_num <- cond_num
	  se_reliable <- !(min_eig <= 0 || cond_num > 1e10)

	  if (min_eig <= 0 || cond_num > 1e10) {
		  cli::cli_warn("Hessian is ill-conditioned (min eigenvalue: {.val {sprintf('%.2e', min_eig)}}, condition number: {.val {sprintf('%.2e', cond_num)}}). Applying ridge regularization. Consider bootstrap SEs via {.fn boot_sir} for reliable inference.")
		  # ridge regularization
		  lambda <- max_eig * 1e-6
		  H_reg <- H + lambda * diag(nrow(H))
		  H_inv <- tryCatch(solve(H_reg), error = function(e) {
			  cli::cli_warn("Regularized Hessian inversion failed. Using {.fn MASS::ginv}.")
			  tryCatch(MASS::ginv(H_reg), error = function(e2) {
				  cli::cli_warn("Generalized inverse also failed. Cannot compute standard errors.")
				  return(NULL)
			  })
		  })
	  } else {
		  H_inv <- tryCatch(solve(H), error = function(e) {
			  cli::cli_warn("Hessian inversion failed unexpectedly. Using {.fn MASS::ginv}.")
			  tryCatch(MASS::ginv(H), error = function(e2) {
				  return(NULL)
			  })
		  })
	  }

	  if (!is.null(H_inv)) {
		  # classical SE
		  se_diag <- diag(H_inv)
		  # ensure non-negative variance
		  se_diag[se_diag < 0] <- NA
		  se <- sqrt(se_diag)

		  # robust (sandwich) SE
		  vcov_mat <- H_inv
		  vcov_robust <- tryCatch(H_inv %*% S %*% H_inv, error = function(e) NULL)

		  rse <- tryCatch({
			  rse_diag <- diag(vcov_robust)
			  rse_diag[rse_diag < 0] <- NA
			  sqrt(rse_diag)
		  }, error = function(e) {
			  cli::cli_warn("Sandwich estimator calculation failed.")
			  return(rep(NA, length(se)))
		  })

		  summ$se <- se
		  summ$rse <- rse
		  summ$t_se <- summ$coef / summ$se
		  summ$t_rse <- summ$coef / summ$rse
	  } else {
		  se_reliable <- FALSE
		  summ$se <- NA
		  summ$rse <- NA
		  summ$t_se <- NA
		  summ$t_rse <- NA
	  }
	}
	}

	# parameter names
	Z_names <- if (q > 0 && !is.null(dimnames(Z)[[3]])) dimnames(Z)[[3]] else paste0("Z", 1:q)
	W_names <- if (p > 0 && !is.null(dimnames(W)[[3]])) dimnames(W)[[3]] else paste0("W", 1:p)

	theta_names <- if (q>0) paste0("(Z) ", Z_names) else NULL

	if (fix_receiver && p > 0) {
	alpha_names <- paste0("(alphaW) ", W_names)
	beta_names <- NULL
	} else {
	alpha_names <- if (p>1) paste0("(alphaW) ", W_names[-1]) else NULL
	beta_names  <- if (p>0) paste0("(betaW) ", W_names) else NULL
	}

	rownames(summ) <- c(theta_names, alpha_names, beta_names)

	# apply names to vcov matrices
	param_names <- rownames(summ)
	if (!is.null(vcov_mat)) {
	  rownames(vcov_mat) <- colnames(vcov_mat) <- param_names
	}
	if (!is.null(vcov_robust)) {
	  rownames(vcov_robust) <- colnames(vcov_robust) <- param_names
	}

	# residuals
	response_resid <- Y - fitted_values

	pearson_resid <- response_resid
	if (family == "poisson") {
	  pearson_resid <- response_resid / sqrt(fitted_values)
	} else if (family == "binomial") {
	  pearson_resid <- response_resid / sqrt(fitted_values * (1 - fitted_values))
	} else if (family == "normal") {
	  pearson_resid <- response_resid / sqrt(sigma2_est)
	}

	deviance_resid <- response_resid
	if (family == "poisson") {
	  dev_contrib <- 2 * (ifelse(Y > 0, Y * log(Y / fitted_values), 0) - (Y - fitted_values))
	  dev_contrib[dev_contrib < 0] <- 0
	  deviance_resid <- sign(response_resid) * sqrt(dev_contrib)
	} else if (family == "binomial") {
	  dev_contrib <- 2 * (ifelse(Y > 0, Y * log(Y / fitted_values), 0) +
						  ifelse(Y < 1, (1 - Y) * log((1 - Y) / (1 - fitted_values)), 0))
	  dev_contrib[dev_contrib < 0] <- 0
	  deviance_resid <- sign(response_resid) * sqrt(dev_contrib)
	}

	resid_list <- list(
	  response = response_resid,
	  pearson = pearson_resid,
	  deviance = deviance_resid
	)

		# count observations (exclude diagonal only for one-mode square networks)
		if (one_mode_square) {
	  N_obs <- 0L
	  for (t in seq_len(T_len)) {
		  Yt <- Y[,,t]
		  diag(Yt) <- NA
		  N_obs <- N_obs + sum(!is.na(Yt))
	  }
	} else {
	  N_obs <- sum(!is.na(Y))
	}
	
	# prepare output
	result <- list(
	summ = summ,
	A    = A,
	B    = B,
	ll   = ll,
	family = family,
	method = method,
	tab = tab,
	theta = theta,
	alpha = alpha,
	beta = beta,
	p = p,
	q = q,
	m = m,
	n1 = n1,
	n2 = n2,
	bipartite = bipartite,
	n_periods = T_len,
	nobs = N_obs,
	fitted.values = fitted_values,
	residuals = resid_list,
	vcov = vcov_mat,
	vcov_robust = vcov_robust,
	Y = Y,
	W = W,
	X = X,
	Z = Z,
	fix_receiver = fix_receiver,
	symmetric = symmetric,
	kron_mode = kron_mode,
	dynamic_W = dynamic_W,
	iterations = if (!is.null(mod$iterations)) mod$iterations else NA,
	history = list(ALPHA=mod$ALPHA, BETA=mod$BETA, THETA=mod$THETA, DEV=mod$DEV),
	convergence = fit_converged,
	se_reliable = se_reliable,
	hessian_cond_num = hessian_cond_num,
	hessian_min_eig = hessian_min_eig,
	seed = seed,
	sigma2_mle = sigma2_mle,
	call = match.call()
	)

	if (family == "normal") {
	  result$sigma2 <- sigma2_est
	}

	# keep package plot dispatch ahead of igraph dispatch
	class(result) <- c("sir_fit", "sir")
	return(result)
}
