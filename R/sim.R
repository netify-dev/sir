
#' Simulate Data from a Social Influence Regression Model
#'
#' Generates synthetic network data from a known SIR data-generating process.
#' Useful for testing, benchmarking, and pedagogical demonstrations. The
#' function simulates Y[i,j,t] from the specified family using influence
#' covariates W, lagged network state X, and optional exogenous covariates Z.
#'
#' @details
#' The influence-carrying state \code{X} is scaled by \code{1 / (m - 1)} before
#' it enters the bilinear mean. The bilinear term \eqn{A X_t B^\top} sums over all
#' \eqn{(m-1)} off-diagonal partners on each side, so without this scaling the
#' linear predictor would grow with network size and (for the Poisson/log link)
#' the conditional mean would saturate, making the influence parameters
#' unrecoverable. The returned \code{X} already includes this scaling, so a plain
#' \code{\link{sir}} fit on \code{(Y, W, X, Z)} recovers the same \code{A}, \code{B}
#' used to generate the data.
#'
#' @param m Integer. Number of nodes in the network.
#' @param T_len Integer. Number of time periods.
#' @param p Integer. Number of influence covariates in W. Default is 2.
#' @param q Integer. Number of exogenous covariates in Z. Default is 1.
#'   Set to 0 for no exogenous covariates.
#' @param family Character string: \code{"poisson"} (default), \code{"normal"},
#'   or \code{"binomial"}.
#' @param alpha Numeric vector of length p for sender influence weights.
#'   The first element (alpha_1) is fixed at 1 for identifiability; only
#'   alpha_2:p are free. If NULL (default), drawn from N(0, 0.3).
#'   Use \code{seed} for reproducibility.
#' @param beta Numeric vector of length p for receiver influence weights.
#'   If NULL (default), drawn from N(0, 0.3).
#' @param theta Numeric vector of length q for exogenous covariate effects.
#'   If NULL (default), drawn from U(-0.5, 0.5).
#' @param W Optional 3D array (m x m x p) of influence covariates. If NULL
#'   (default), generated with standard normal entries (a dense, well-conditioned
#'   influence design). The diagonal is set to zero because self-ties are not
#'   part of the one-mode SIR likelihood.
#' @param sigma Numeric. Standard deviation for the normal family. Default 1.
#' @param seed Optional integer for reproducibility. When supplied, the seed is
#'   set locally and the caller's global RNG state is restored on exit, so a
#'   subsequent draw (e.g. \code{runif}) in the caller is left unperturbed.
#' @param ... Unused; catches mistyped arguments and reports a clear error.
#' @return A list with components:
#'   \describe{
#'     \item{Y}{3D array (m x m x T_len) of simulated outcomes.}
#'     \item{W}{3D array (m x m x p) of influence covariates.}
#'     \item{X}{3D array (m x m x T_len) of the (scaled) lagged network state
#'       used in the bilinear mean: \code{X[,,t]} is \code{log(Y[,,t-1] + 1)}
#'       (Poisson) or \code{Y[,,t-1]} (otherwise), divided by \code{(m - 1)}.}
#'     \item{Z}{4D array (m x m x q x T_len) of exogenous covariates, or
#'       NULL if q = 0.}
#'     \item{alpha}{True alpha vector (length p, with alpha_1 = 1).}
#'     \item{beta}{True beta vector (length p).}
#'     \item{theta}{True theta vector (length q).}
#'     \item{A}{True sender influence matrix (m x m).}
#'     \item{B}{True receiver influence matrix (m x m).}
#'     \item{family}{The distribution family used.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Simulate Poisson network and recover the parameters
#' dat <- sim_sir(m = 15, T_len = 30, p = 2, q = 1, family = "poisson", seed = 42)
#' fit <- sir(dat$Y, dat$W, dat$X, dat$Z, family = "poisson")
#' cbind(true = c(dat$theta, dat$alpha[-1], dat$beta), estimated = coef(fit))
#' }
#' @export
sim_sir <- function(m, T_len, p = 2, q = 1, family = "poisson",
					alpha = NULL, beta = NULL, theta = NULL,
					W = NULL, sigma = 1, seed = NULL, ...) {

	# catch mistyped arguments with a helpful message (e.g. time= or T= for T_len)
	if (...length() > 0) {
		bad_args <- ...names()
		cli::cli_abort(c(
			"Unknown argument{?s} to {.fn sim_sir}: {.val {bad_args}}.",
			"i" = "The number of time periods is {.arg T_len}; see {.code ?sim_sir} for all arguments."
		))
	}

	# reproducible draws without leaking state: when a seed is supplied we set it
	# locally and restore the caller's RNG stream on exit, so a later runif() in
	# the caller is unaffected (mirrors sir()).
	if (!is.null(seed)) {
		if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
			old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
			on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
		} else {
			on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
				rm(".Random.seed", envir = globalenv()), add = TRUE)
		}
		set.seed(seed)
	}

	if (!family %in% c("poisson", "normal", "binomial")) {
		cli::cli_abort("family must be one of {.val poisson}, {.val normal}, {.val binomial}.")
	}

	# generate or validate parameters
	if (is.null(alpha)) {
		alpha <- c(1, rnorm(max(p - 1, 0), sd = 0.3))
	} else {
		if (length(alpha) != p) cli::cli_abort("{.arg alpha} must have length {.val {p}}.")
		alpha[1] <- 1
	}

	beta_was_null <- is.null(beta)
	if (is.null(beta)) {
		beta <- rnorm(p, sd = 0.3)
	} else {
		if (length(beta) != p) cli::cli_abort("{.arg beta} must have length {.val {p}}.")
	}

	if (q > 0) {
		if (is.null(theta)) {
			theta <- runif(q, -0.5, 0.5)
		} else {
			if (length(theta) != q) cli::cli_abort("{.arg theta} must have length {.val {q}}.")
		}
	} else {
		theta <- numeric(0)
	}

	# generate dense influence covariates
	if (is.null(W)) {
		W <- array(rnorm(m * m * p), dim = c(m, m, p))
	} else {
		if (!all(dim(W)[1:2] == m) || dim(W)[3] != p) {
			cli::cli_abort("W dimensions must be {.val {m}} x {.val {m}} x {.val {p}}.")
		}
	}
	W <- set_square_diagonal(W, 0)

	# build influence matrices
	A <- matrix(0, m, m)
	B <- matrix(0, m, m)
	if (p > 0) {
		for (k in seq_len(p)) {
			A <- A + alpha[k] * W[, , k]
			B <- B + beta[k] * W[, , k]
		}
	}

	# scale the influence-carrying state so the linear predictor stays O(1) across
	# network sizes; baked into X so a plain sir() fit recovers the same A, B
	infl_scale <- max(m - 1, 1)

	# keep the count/continuous recursion stationary by rescaling beta when the
	# bilinear gain exceeds 1; only auto-generated beta is touched, binomial is bounded
	if (beta_was_null && family %in% c("poisson", "normal")) {
		gain <- max(abs(eigen(A, only.values = TRUE)$values)) *
				max(abs(eigen(B, only.values = TRUE)$values)) / infl_scale
		if (is.finite(gain) && gain > 0.8) {
			scale_b <- 0.8 / gain
			beta <- beta * scale_b
			B <- B * scale_b
		}
	}

	# generate Z
	Z <- NULL
	if (q > 0) {
		Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	}

	# generate Y sequentially (X_t carries the lagged state Y_{t-1})
	Y <- array(0, dim = c(m, m, T_len))
	X <- array(0, dim = c(m, m, T_len))

	# initialize Y[,,1] from an intercept-only draw
	if (family == "poisson") {
		Y[,,1] <- matrix(rpois(m * m, lambda = 2), m, m)
	} else if (family == "normal") {
		Y[,,1] <- matrix(rnorm(m * m, sd = sigma), m, m)
	} else {
		Y[,,1] <- matrix(rbinom(m * m, 1, 0.3), m, m)
	}
	diag(Y[,,1]) <- 0

	for (t in seq_len(T_len)[-1]) {
		# influence-carrying state from the lag; log-transform counts so the
		# exp() link does not produce explosive dynamics. then apply the scaling.
		if (family == "poisson") {
			x_raw <- log(Y[,,t - 1] + 1)
		} else {
			x_raw <- Y[,,t - 1]
		}
		x_raw[is.na(x_raw)] <- 0
		X[,,t] <- x_raw / infl_scale

		# build the linear predictor
		eta <- A %*% X[,,t] %*% t(B)
		if (q > 0) {
			for (k in 1:q) {
				eta <- eta + theta[k] * Z[,,k,t]
			}
		}

		# generate from the family
		if (family == "poisson") {
			lambda <- exp(eta)
			# generous guard against numerical overflow; with the scaling above
			# this effectively never binds for reasonable coefficients.
			lambda[lambda > 1e6] <- 1e6
			lambda[lambda < 1e-10] <- 1e-10
			Y[,,t] <- matrix(rpois(m * m, lambda = c(lambda)), m, m)
		} else if (family == "normal") {
			Y[,,t] <- matrix(rnorm(m * m, mean = c(eta), sd = sigma), m, m)
		} else {
			prob <- 1 / (1 + exp(-eta))
			Y[,,t] <- matrix(rbinom(m * m, 1, c(prob)), m, m)
		}
		diag(Y[,,t]) <- 0
	}

	list(
		Y = Y, W = W, X = X, Z = Z,
		alpha = alpha, beta = beta, theta = theta,
		A = A, B = B, family = family
	)
}
