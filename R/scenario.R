#' Get Scenario Values for Prediction
#'
#' Computes representative covariate values for model-implied scenario analysis.
#' For each influence covariate, returns a set of values (the 10th,
#' 25th, 50th, 75th, and 90th percentiles) at which to evaluate predicted
#' influence.
#'
#' @param data A 3D (\code{m x m x p}) or 4D (\code{m x m x p x T}) array of
#'   influence covariates, or a data frame of covariates.
#' @param vars Character vector of variable names to compute scenario values
#'   for. If NULL, uses all available variables.
#' @param time Optional time index for time-varying covariates. If NULL, values
#'   are pooled across time.
#' @param directed Logical. If TRUE (default), treats the network as directed
#'   when computing scenario values. For square arrays, diagonals are excluded
#'   unless \code{include_diag = TRUE}; if FALSE, only the upper triangle is used.
#' @param include_diag Logical. If TRUE, include square-array diagonal cells when
#'   computing quantiles. Keep the default FALSE for one-mode networks; use TRUE
#'   for square bipartite sender/receiver arrays whose diagonal cells are real.
#' @param node_names Optional character vector of node names for labeling.
#' @param scen_vals Optional named list of pre-computed scenario values.
#'
#' @return A named list mapping each variable to a vector of quantile values.
#'
#' @examples
#' dat <- sim_sir(m = 8, T_len = 10, p = 2, q = 1, family = "poisson", seed = 1)
#' sv <- get_scen_vals(dat$W)
#' str(sv)
#' @export
get_scen_vals <- function(
	data,
	vars = NULL,
		time = NULL,
		directed = TRUE,
		include_diag = FALSE,
		node_names = NULL,
		scen_vals = NULL
) {
	if (!is.null(scen_vals)) {
		if (!is.null(vars)) scen_vals <- scen_vals[vars]
		return(scen_vals)
	}

	qfun <- function(vals) {
		vals <- vals[is.finite(vals)]
		if (!length(vals)) {
			cli::cli_abort("No finite values available for scenario quantiles.")
		}
		stats::quantile(vals, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
	}

	if (is.array(data)) {
		d <- dim(data)
		if (!(length(d) %in% c(3L, 4L))) {
			cli::cli_abort("{.arg data} must be a 3D or 4D array, or a data frame.")
		}
		var_names <- dimnames(data)[[3]]
		if (is.null(var_names)) var_names <- paste0("var", seq_len(d[3]))
		use_vars <- if (is.null(vars)) var_names else vars
		missing_vars <- setdiff(use_vars, var_names)
		if (length(missing_vars)) {
			cli::cli_abort("{.arg vars} contains unknown variable{?s}: {.val {missing_vars}}.")
		}
		var_idx <- match(use_vars, var_names)
		if (!is.null(time)) {
			if (length(d) != 4L) {
				cli::cli_abort("{.arg time} can only be used with 4D time-varying arrays.")
			}
			if (!is.numeric(time) || !length(time) || any(!is.finite(time)) ||
				any(time < 1) || any(time > d[4]) || any(time != round(time))) {
				cli::cli_abort("{.arg time} must index the fourth dimension of {.arg data}.")
			}
			time <- as.integer(time)
		}
		if (!is.null(node_names) && length(node_names) != d[1]) {
			cli::cli_abort("{.arg node_names} has length {length(node_names)} but {.arg data} has {d[1]} rows.")
		}

			cell_mask <- matrix(TRUE, d[1], d[2])
			if (d[1] == d[2]) {
				if (isTRUE(directed)) {
					if (!isTRUE(include_diag)) {
						cell_mask[row(cell_mask) == col(cell_mask)] <- FALSE
					}
				} else {
					cell_mask <- upper.tri(cell_mask, diag = isTRUE(include_diag))
				}
			}

		scenario <- list()
		for (ii in seq_along(var_idx)) {
			vi <- var_idx[ii]
			vals <- if (length(d) == 4L) {
				t_idx <- if (is.null(time)) seq_len(d[4]) else time
				unlist(lapply(t_idx, function(tt) data[, , vi, tt][cell_mask]), use.names = FALSE)
			} else {
				data[, , vi][cell_mask]
			}
			scenario[[use_vars[ii]]] <- qfun(as.numeric(vals))
		}
		return(scenario)
	}

	if (is.data.frame(data)) {
		use_vars <- if (is.null(vars)) names(data) else vars
		missing_vars <- setdiff(use_vars, names(data))
		if (length(missing_vars)) {
			cli::cli_abort("{.arg vars} contains unknown variable{?s}: {.val {missing_vars}}.")
		}
		scenario <- list()
		for (nm in use_vars) {
			if (!is.numeric(data[[nm]])) {
				cli::cli_abort("{.arg data} column {.val {nm}} is not numeric.")
			}
			scenario[[nm]] <- qfun(as.numeric(data[[nm]]))
		}
		return(scenario)
	}

	if (!is.null(vars)) {
		cli::cli_abort("{.arg data} must be an array or data frame unless {.arg scen_vals} is supplied.")
	}

	list()
}

#' Build Scenario Array for Prediction
#'
#' Constructs an artificial influence covariate grid for model-implied scenario prediction,
#' where one variable is set to a scenario value while the others are held at the
#' mean of their supplied scenario values. The resulting array can be passed to
#' \code{\link{predict.sir}} via \code{newdata}. The diagonal is set to zero by
#' default for one-mode scenario arrays so self-influence channels are not
#' reintroduced. For empirical counterfactuals, prefer copying the observed
#' \code{W} array and modifying a theoretically defined subset of cells.
#'
#' @param var_to_vary Character name of the variable to vary.
#' @param scen_vals Named list of scenario values (from \code{\link{get_scen_vals}}).
#' @param node_names Character vector of node names.
#' @param var_names Character vector of all variable names in the array.
#' @param n_time Optional integer number of time periods. If NULL (default),
#'   returns a static 3D \code{n x n x p} array. Set this to return a 4D
#'   time-varying \code{n x n x p x n_time} array.
#' @param value Numeric scenario value for \code{var_to_vary}. If NULL, uses the
#'   mean of that variable's scenario values.
#' @param zero_diag Logical. If TRUE (default), set square-array diagonals to
#'   zero. Use FALSE for square bipartite sender-side arrays whose diagonal cells
#'   are real observations.
#'
#' @return A static 3D scenario array by default, or a 4D array when
#'   \code{n_time} is supplied.
#'
#' @examples
#' dat <- sim_sir(m = 8, T_len = 10, p = 2, q = 1, family = "poisson", seed = 1)
#' sv <- get_scen_vals(dat$W)
#' arr <- get_scen_array(names(sv)[1], scen_vals = sv,
#'                       node_names = paste0("n", 1:8), var_names = names(sv))
#' dim(arr)
#' @export
get_scen_array <- function(
	var_to_vary,
	scen_vals,
	node_names,
	var_names,
	n_time = NULL,
	value = NULL,
	zero_diag = TRUE
) {
	n <- length(node_names)
	p <- length(var_names)
	if (!is.null(n_time) &&
		(!is.numeric(n_time) || length(n_time) != 1L || n_time < 1 ||
			!is.finite(n_time) || n_time != round(n_time))) {
		cli::cli_abort("{.arg n_time} must be NULL or a positive integer.")
	}
	if (is.null(n_time)) {
		scen <- array(
			0,
			dim = c(n, n, p),
			dimnames = list(node_names, node_names, var_names)
		)
	} else {
		n_time <- as.integer(n_time)
		scen <- array(
			0,
			dim = c(n, n, p, n_time),
			dimnames = list(node_names, node_names, var_names, NULL)
		)
	}

	idx <- match(var_to_vary, var_names)
	if (is.na(idx)) {
		cli::cli_abort("{.arg var_to_vary} not found in {.arg var_names}.")
	}
	if (is.null(scen_vals[[var_to_vary]])) {
		cli::cli_abort("{.arg scen_vals} has no values for {.val {var_to_vary}}.")
	}
	target_value <- if (is.null(value)) {
		mean(scen_vals[[var_to_vary]], na.rm = TRUE)
	} else {
		value
	}
	if (!is.numeric(target_value) || length(target_value) != 1 || !is.finite(target_value)) {
		cli::cli_abort("{.arg value} must be a single finite number.")
	}

	# fill reference values
	for (vi in seq_len(p)) {
		ref_value <- if (!is.null(scen_vals[[var_names[vi]]])) {
			mean(scen_vals[[var_names[vi]]])
		} else {
			0
		}
		if (length(dim(scen)) == 3L) {
			scen[, , vi] <- ref_value
		} else {
			scen[, , vi, ] <- ref_value
		}
	}

	# set the target scenario value
	if (length(dim(scen)) == 3L) {
		scen[, , idx] <- target_value
	} else {
		scen[, , idx, ] <- target_value
	}
	if (isTRUE(zero_diag)) scen <- set_square_diagonal(scen, 0)
	scen
}
