
#' Construct Relational Covariates from a Network Array
#'
#' Builds standard relational covariates from a base network array: the
#' original (main) effect, the reciprocal (transpose) effect, and a
#' transitive closure effect. These are common exogenous covariates (Z)
#' in SIR models capturing higher-order network dependencies.
#'
#' @param arr A 3D array (m x m x T) of network data. Typically an outcome
#'   variable or a covariate that varies across dyads and time.
#' @param name Character string used to label the covariates in the output
#'   array's third dimension. The main effect is labeled \code{name}, the
#'   reciprocal \code{paste0(name, "_recip")}, and the transitive
#'   \code{paste0(name, "_trans")}.
#' @param effects Character vector specifying which relational covariates
#'   to include. Any subset of \code{c("main", "reciprocal", "transitive")}.
#'   Default is all three.
#'
#' @details
#' \describe{
#'   \item{Main}{The original array, Z_ij = arr_ij. This captures the
#'     direct dyadic effect.}
#'   \item{Reciprocal}{The transpose, Z_ij = arr_ji. Captures whether the
#'     reverse relationship matters.}
#'   \item{Transitive}{A measure of shared connectivity:
#'     Z_ij = (S %*% S)_ij where S = (arr + arr')/2 is the symmetrized
#'     network. Captures triadic closure and transitivity effects.}
#' }
#'
#' @return A 4D array (m x m x q x T) where q is the number of requested
#'   effects (1 to 3). Suitable for passing directly as the \code{Z}
#'   argument to \code{\link{sir}}.
#'
#' @examples
#' \dontrun{
#' # Build relational covariates from trade data
#' Z_trade <- rel_covar(trade_array, "trade")
#' dim(Z_trade)  # m x m x 3 x T
#'
#' # Use only main and reciprocal effects
#' Z_simple <- rel_covar(trade_array, "trade", effects = c("main", "reciprocal"))
#'
#' # Pass to sir() as exogenous covariates
#' fit <- sir(Y, W, X, Z = Z_trade, family = "poisson")
#' }
#' @export
rel_covar <- function(arr, name, effects = c("main", "reciprocal", "transitive")) {

	effects <- match.arg(effects, c("main", "reciprocal", "transitive"),
						 several.ok = TRUE)

	dims <- dim(arr)
	if (length(dims) != 3) {
		cli::cli_abort("{.arg arr} must be a 3D array (m x m x T), got {.val {length(dims)}}D.")
	}

	m <- dims[1]
	T_len <- dims[3]

	# build each requested effect
	covs <- list()
	cov_names <- character(0)

	if ("main" %in% effects) {
		covs[[length(covs) + 1]] <- arr
		cov_names <- c(cov_names, name)
	}

	if ("reciprocal" %in% effects) {
		covs[[length(covs) + 1]] <- aperm(arr, c(2, 1, 3))
		cov_names <- c(cov_names, paste0(name, "_recip"))
	}

	if ("transitive" %in% effects) {
		trans <- array(0, dim = dims)
		for (t in seq_len(T_len)) {
			S <- (arr[,,t] + t(arr[,,t])) / 2
			trans[,,t] <- S %*% S
		}
		covs[[length(covs) + 1]] <- trans
		cov_names <- c(cov_names, paste0(name, "_trans"))
	}

	q <- length(covs)
	result <- array(NA, dim = c(m, m, q, T_len))

	for (k in seq_len(q)) {
		result[,,k,] <- covs[[k]]
	}

	# set dimension names
	dn <- dimnames(arr)
	dimnames(result) <- list(dn[[1]], dn[[2]], cov_names, dn[[3]])

	result
}
