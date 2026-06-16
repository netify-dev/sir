
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
#' For the Poisson and Normal families the lagged recursion is stationary only
#' when the spectral gain \eqn{\rho(A)\rho(B)/(m-1)} is below 1. By default,
#' auto-generated coefficients are rescaled to a conservative gain (at most 0.8) so
#' simulations are reliably stable, and user-supplied coefficients are used exactly
#' as given (a gain \eqn{\ge 1} triggers an explosive-series warning). Set
#' \code{gain} to target a specific value: the influence operators are then
#' rescaled so \eqn{\rho(A)\rho(B)/(m-1)} equals \code{gain} exactly, letting the
#' bilinear term carry a chosen, strong-but-stationary share of the dynamics. A
#' larger \code{gain} (say 0.9) makes the influence mechanism dominate the lagged
#' dynamics while still recovering cleanly; values near 1 approach
#' non-stationarity. Binomial outcomes are bounded, so \code{gain} does not apply.
#'
#' @param m Integer. Number of nodes in the network.
#' @param T_len Integer. Number of time periods.
#' @param p Integer. Number of influence covariates in W. Default is 2.
#' @param q Integer. Number of exogenous covariates in Z. Default is 1.
#'   Set to 0 for no exogenous covariates.
#' @param family Character string: \code{"poisson"} (default), \code{"normal"},
#'   or \code{"binomial"}.
#' @param alpha Numeric vector of length p for sender influence weights. For a
#'   directed process the first element (alpha_1) is fixed at 1 for
#'   identifiability; for \code{symmetric = TRUE} all elements are kept as given
#'   (the anchor need not be 1). If NULL (default), drawn from N(0, 0.3) with
#'   alpha_1 = 1. Use \code{seed} for reproducibility.
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
#' @param symmetric Logical. If TRUE, simulate a genuine \strong{undirected}
#'   network: \code{W} is symmetrized, \code{beta} is tied to \code{alpha} so
#'   \code{B = A}, and each \code{Y}/\code{Z} slice is symmetric. For the
#'   count/continuous recursion the influence covariates are rescaled (when
#'   auto-generated) so the symmetric gain \eqn{\rho(A)^2/(m-1)} stays below 1.
#'   Pairs with \code{sir(..., symmetric = TRUE)}, which fits and reports the
#'   shared operator as \code{gamma}; \code{sim_sir} stores that same vector in
#'   \code{$alpha} (which equals \code{$beta} here). Default FALSE.
#' @param gain Optional numeric in (0, 1). Target spectral gain
#'   \eqn{\rho(A)\rho(B)/(m-1)} for the Poisson/Normal lagged recursion. When
#'   supplied, the influence operators are rescaled to hit it exactly (scaling
#'   \code{beta}/\code{B} for directed fits, the shared operator for symmetric,
#'   leaving \code{alpha_1 = 1} intact), so the influence term carries a chosen,
#'   strong-but-stationary share of the dynamics. Use a larger value (e.g. 0.9)
#'   for a strong, recoverable influence signal; values near 1 approach
#'   non-stationarity. NULL (default) keeps the conservative auto-rescaling.
#'   Ignored for \code{family = "binomial"}.
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
#'     \item{alpha}{True alpha vector (length p; alpha_1 = 1 for directed, kept
#'       as supplied for symmetric). For symmetric fits this equals \code{beta}
#'       and is the shared operator reported by the fit as gamma.}
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
					W = NULL, sigma = 1, seed = NULL, symmetric = FALSE,
					gain = NULL, ...) {

	# target spectral gain: rescale the influence operators so the bilinear term
	# carries a chosen, strong-but-stationary share of the dynamics
	if (!is.null(gain)) {
		if (!is.numeric(gain) || length(gain) != 1 || !is.finite(gain) ||
			gain <= 0 || gain >= 1) {
			cli::cli_abort(c(
				"{.arg gain} must be a single number in (0, 1).",
				"i" = "It sets the target spectral gain rho(A)rho(B)/(m-1); values near 1 give strong but near-non-stationary influence."
			))
		}
		if (family == "binomial") {
			cli::cli_warn("{.arg gain} has no effect for the {.val binomial} family (bounded outcomes); ignoring it.")
		}
	}

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

	# generate or validate parameters. directed fits fix alpha_1 = 1 for
	# identifiability; symmetric fits estimate all gamma, so an explicit non-unit
	# anchor is kept as given
	alpha_was_explicit <- !is.null(alpha)
	if (is.null(alpha)) {
		alpha <- c(1, rnorm(max(p - 1, 0), sd = 0.3))
	} else {
		if (length(alpha) != p) cli::cli_abort("{.arg alpha} must have length {.val {p}}.")
		if (!isTRUE(symmetric)) alpha[1] <- 1
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
	W_was_null <- is.null(W)
	if (is.null(W)) {
		W <- array(rnorm(m * m * p), dim = c(m, m, p))
	} else {
		if (!all(dim(W)[1:2] == m) || dim(W)[3] != p) {
			cli::cli_abort("W dimensions must be {.val {m}} x {.val {m}} x {.val {p}}.")
		}
	}
	W <- set_square_diagonal(W, 0)

	# undirected: symmetrize W and tie beta = alpha so B = A
	if (isTRUE(symmetric)) {
		for (k in seq_len(p)) W[, , k] <- 0.5 * (W[, , k] + t(W[, , k]))
		W <- set_square_diagonal(W, 0)
		beta <- alpha
	}

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

	# keep the count/continuous recursion stationary. auto-generated coefficients
	# are rescaled so the bilinear gain stays below 1; user-supplied coefficients
	# are left exactly as given, but a gain >= 1 (an explosive process) is warned.
	# binomial is bounded, so the gain does not apply.
	if (family %in% c("poisson", "normal")) {
		rho_A <- max(abs(eigen(A, only.values = TRUE)$values))
		rho_B <- if (isTRUE(symmetric)) rho_A else max(abs(eigen(B, only.values = TRUE)$values))
		gain_cur <- rho_A * rho_B / infl_scale
		auto <- if (isTRUE(symmetric)) W_was_null else beta_was_null
		if (!is.null(gain) && is.finite(gain_cur) && gain_cur > 0) {
			# rescale to hit the requested gain exactly (scale beta/B for directed,
			# the shared operator for symmetric, leaving alpha_1 = 1 intact)
			if (isTRUE(symmetric)) {
				s <- sqrt(gain / gain_cur); W <- W * s; A <- A * s; B <- A
			} else {
				s <- gain / gain_cur; beta <- beta * s; B <- B * s
			}
		} else if (is.null(gain) && is.finite(gain_cur) && gain_cur > 0.8 && auto) {
			if (isTRUE(symmetric)) {
				s <- sqrt(0.8 / gain_cur); W <- W * s; A <- A * s; B <- A
			} else {
				scale_b <- 0.8 / gain_cur; beta <- beta * scale_b; B <- B * scale_b
			}
		} else if (is.null(gain) && is.finite(gain_cur) && gain_cur >= 1 && !auto) {
			cli::cli_warn(c(
				"The supplied influence coefficients imply a non-stationary process (spectral gain {.val {sprintf('%.2f', gain_cur)}} >= 1).",
				"i" = "The simulated {.arg Y} can grow explosively; use smaller coefficients or set {.arg gain} for a stationary series."
			))
		}
	}

	# generate Z (symmetrize each slice for the undirected DGP)
	Z <- NULL
	if (q > 0) {
		Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
		if (isTRUE(symmetric)) {
			for (k in seq_len(q)) for (t in seq_len(T_len)) {
				Z[, , k, t] <- 0.5 * (Z[, , k, t] + t(Z[, , k, t]))
			}
		}
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
	if (isTRUE(symmetric)) {
		Y1 <- Y[,,1]
		if (family == "normal") Y1 <- 0.5 * (Y1 + t(Y1)) else Y1[lower.tri(Y1)] <- t(Y1)[lower.tri(Y1)]
		Y[,,1] <- Y1
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
		# symmetrize the draw; copy the upper triangle down for discrete families
		# so values stay valid, average for normal
		if (isTRUE(symmetric)) {
			Yt <- Y[,,t]
			if (family == "normal") {
				Yt <- 0.5 * (Yt + t(Yt))
			} else {
				# copy upper triangle to lower so counts/binary stay valid
				Yt[lower.tri(Yt)] <- t(Yt)[lower.tri(Yt)]
			}
			Y[,,t] <- Yt
		}
		diag(Y[,,t]) <- 0
	}

	list(
		Y = Y, W = W, X = X, Z = Z,
		alpha = alpha, beta = beta, theta = theta,
		A = A, B = B, family = family, symmetric = isTRUE(symmetric)
	)
}
