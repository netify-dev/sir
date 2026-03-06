
#' Simulate Data from a Social Influence Regression Model
#'
#' Generates synthetic network data from a known SIR data-generating process.
#' Useful for testing, benchmarking, and pedagogical demonstrations. The
#' function simulates Y[i,j,t] from the specified family using influence
#' covariates W, lagged network state X, and optional exogenous covariates Z.
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
#' @param beta Numeric vector of length p for receiver influence weights.
#'   If NULL (default), drawn from N(0, 0.3).
#' @param theta Numeric vector of length q for exogenous covariate effects.
#'   If NULL (default), drawn from U(-0.5, 0.5).
#' @param W Optional 3D array (m x m x p) of influence covariates. If NULL
#'   (default), generated with standard normal entries.
#' @param sigma Numeric. Standard deviation for the normal family. Default 1.
#' @param seed Optional integer for reproducibility.
#' @return A list with components:
#'   \describe{
#'     \item{Y}{3D array (m x m x T_len) of simulated outcomes.}
#'     \item{W}{3D array (m x m x p) of influence covariates.}
#'     \item{X}{3D array (m x m x T_len) of lagged network state.}
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
#' # Simulate Poisson network
#' dat <- sim_sir(m = 15, T_len = 10, p = 2, q = 1, family = "poisson", seed = 42)
#'
#' # Fit model to recover parameters
#' fit <- sir(dat$Y, dat$W, dat$X, dat$Z, family = "poisson")
#' cbind(true = c(dat$theta, dat$alpha[-1], dat$beta), estimated = coef(fit))
#' }
#' @export
sim_sir <- function(m, T_len, p = 2, q = 1, family = "poisson",
					alpha = NULL, beta = NULL, theta = NULL,
					W = NULL, sigma = 1, seed = NULL) {

	if (!is.null(seed)) set.seed(seed)

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

	# generate W
	if (is.null(W)) {
		W <- array(rnorm(m * m * p), dim = c(m, m, p))
		# make W[,,1] = identity-like (common in practice)
		W[,,1] <- diag(m) + W[,,1] * 0.1
	} else {
		if (!all(dim(W)[1:2] == m) || dim(W)[3] != p) {
			cli::cli_abort("W dimensions must be {.val {m}} x {.val {m}} x {.val {p}}.")
		}
	}

	# construct A, B
	A <- matrix(0, m, m)
	B <- matrix(0, m, m)
	for (k in 1:p) {
		A <- A + alpha[k] * W[,,k]
		B <- B + beta[k] * W[,,k]
	}
	diag(A) <- 0
	diag(B) <- 0

	# generate Z
	Z <- NULL
	if (q > 0) {
		Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	}

	# generate Y sequentially (Y_t depends on X_t = Y_{t-1})
	Y <- array(0, dim = c(m, m, T_len))
	X <- array(0, dim = c(m, m, T_len))

	# initialize Y[,,1] from intercept-only model
	if (family == "poisson") {
		Y[,,1] <- matrix(rpois(m * m, lambda = 2), m, m)
	} else if (family == "normal") {
		Y[,,1] <- matrix(rnorm(m * m, sd = sigma), m, m)
	} else {
		Y[,,1] <- matrix(rbinom(m * m, 1, 0.3), m, m)
	}
	diag(Y[,,1]) <- 0

	for (t in 2:T_len) {
		X[,,t] <- Y[,,t - 1]

		# linear predictor: eta = theta'Z + A X B'
		eta <- A %*% X[,,t] %*% t(B)
		if (q > 0) {
			for (k in 1:q) {
				eta <- eta + theta[k] * Z[,,k,t]
			}
		}

		# generate from family
		if (family == "poisson") {
			lambda <- exp(eta)
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
