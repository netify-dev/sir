
#' Compute Summary Statistics for Scenario Construction
#'
#' Extracts summary statistics (mean, quantiles) from a 4D covariate
#' array.  These values are used to set baseline and varying levels
#' when building counterfactual scenario arrays with
#' \code{\link{get_scen_array}}.
#'
#' @param data Four-dimensional array (m x m x p x T) of covariates,
#'   or a three-dimensional array (m x m x T) treated as a single
#'   variable.
#' @param vars Character vector of variable names (matching the third
#'   dimension of \code{data}) to summarise.  If NULL, all variables
#'   are used.
#' @param time Integer vector of time indices to include.  If NULL, all
#'   time periods are used.
#' @param directed Logical vector of length \code{length(vars)}
#'   indicating whether each variable is directed.
#'   For undirected variables the lower triangle is excluded before
#'   computing statistics.  Defaults to TRUE for all variables.
#'
#' @return A named list.  Each element corresponds to one variable
#'   and contains:
#'   \describe{
#'     \item{mean}{Scalar mean of off-diagonal entries.}
#'     \item{quantiles}{Named numeric vector of quantiles at 5\%
#'       increments from 0 to 1.}
#'   }
#'
#' @examples
#' \dontrun{
#' # W is m x m x p (or m x m x p x T)
#' vals <- get_scen_vals(W)
#' vals[["proximity"]]$mean
#' vals[["proximity"]]$quantiles
#' }
#'
#' @export
get_scen_vals <- function(data, vars = NULL, time = NULL,
						  directed = NULL) {
	dims <- dim(data)
	ndim <- length(dims)

	if (ndim == 3) {
		# treat (m x m x p) as p variables at a single time point
		dn <- dimnames(data)
		data <- array(data, dim = c(dims[1], dims[2], dims[3], 1),
					  dimnames = list(dn[[1]], dn[[2]], dn[[3]], NULL))
		dims <- dim(data)
	}

	if (length(dims) != 4) {
		cli::cli_abort("{.arg data} must be a 3D or 4D array.")
	}

	m <- dims[1]
	p <- dims[3]
	T_len <- dims[4]

	var_names <- dimnames(data)[[3]]
	if (is.null(var_names)) var_names <- paste0("V", seq_len(p))
	if (is.null(vars)) vars <- var_names
	if (is.null(time)) time <- seq_len(T_len)
	if (is.null(directed)) directed <- rep(TRUE, length(vars))

	vals <- vector("list", length(vars))
	names(vals) <- vars

	for (vi in seq_along(vars)) {
		v <- vars[vi]
		v_idx <- match(v, var_names)
		if (is.na(v_idx)) cli::cli_abort("Variable {.val {v}} not found in data.")

		# collect off-diagonal values
		pieces <- lapply(time, function(t) {
			x <- data[, , v_idx, t]
			diag(x) <- NA
			if (!directed[vi]) x[lower.tri(x)] <- NA
			c(x)
		})
		all_vals <- unlist(pieces)
		all_vals <- all_vals[!is.na(all_vals)]

		vals[[v]] <- list(
			mean = mean(all_vals),
			quantiles = quantile(all_vals, probs = seq(0, 1, 0.05))
		)
	}

	vals
}


#' Build Counterfactual Scenario Array
#'
#' Constructs a 4D array of covariate values for counterfactual
#' prediction, varying one variable across its empirical range while
#' holding others at their mean.
#'
#' @param var_to_vary Character string naming the variable to vary.
#' @param scen_vals Named list produced by \code{\link{get_scen_vals}}.
#' @param node_names Character vector of node names (used for the
#'   first two array dimensions).
#' @param var_names Character vector of all variable names (third
#'   dimension of the output).
#'
#' @return A four-dimensional array of dimension
#'   \code{m x m x p x S}, where \code{S} is the number of
#'   unique quantile values for the varying variable.  Each slice
#'   along the fourth dimension holds one scenario.  Off-diagonal
#'   entries for the varying variable take its quantile value, while
#'   all other variables are set to their mean.  Diagonals are zero.
#'
#' @examples
#' \dontrun{
#' vals <- get_scen_vals(W)
#' scen <- get_scen_array("proximity", vals,
#'                        node_names = rownames(W),
#'                        var_names  = dimnames(W)[[3]])
#' }
#'
#' @export
get_scen_array <- function(var_to_vary, scen_vals, node_names,
						   var_names) {
	m <- length(node_names)
	p <- length(var_names)

	# baseline values: mean for all variables
	base <- vapply(scen_vals[var_names], function(x) x$mean, numeric(1))

	# varying levels: unique quantiles
	vary_vals <- unique(scen_vals[[var_to_vary]]$quantiles)
	S <- length(vary_vals)

	scen <- array(0, dim = c(m, m, p, S),
				  dimnames = list(node_names, node_names, var_names, NULL))

	pos <- match(var_to_vary, var_names)

	for (s in seq_len(S)) {
		vals_s <- base
		vals_s[pos] <- vary_vals[s]
		for (k in seq_len(p)) {
			scen[, , k, s] <- vals_s[k]
		}
		# zero diagonals
		for (k in seq_len(p)) {
			diag(scen[, , k, s]) <- 0
		}
	}

	scen
}
