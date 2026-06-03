# rolling-origin cross-validation and out-of-sample scoring for sir models.
# cv_sir refits sir() on expanding windows and scores the next period via
# forecast(); score_sir computes family-appropriate proper scores.

#' Score predicted networks against observed outcomes
#'
#' Computes out-of-sample scores comparing predicted to actual networks. By
#' default, square arrays are scored on the off-diagonal non-missing cells because
#' one-mode self-ties are not modeled; set \code{drop_diagonal = FALSE} for square
#' bipartite arrays whose diagonal cells are real sender-receiver observations.
#'
#' @param actual,predicted numeric arrays or matrices of identical shape, on the
#'   response scale (counts, probabilities, or means).
#' @param family character; one of \code{"poisson"}, \code{"normal"},
#'   \code{"binomial"}.
#' @param drop_diagonal logical; if TRUE (default), exclude diagonal cells for
#'   square matrices/arrays before scoring.
#' @return a named numeric vector of scores. \code{normal}: rmse, mae.
#'   \code{poisson}: rmse, mae, deviance. \code{binomial}: rmse, mae, logloss,
#'   brier, auc. Lower is better for every score except \code{auc} (higher is
#'   better). Notes: \code{deviance} is the per-cell \emph{mean} unit deviance
#'   (not the total, so it is comparable across origins of different size);
#'   non-finite \code{actual} cells are excluded, while non-finite
#'   \code{predicted} values on scored cells raise an error;
#'   \code{logloss} clamps
#'   predicted probabilities away from 0/1, so against a hard-label (0/1)
#'   baseline it is dominated by that clamp -- prefer \code{brier}/\code{auc}
#'   there; \code{auc} is \code{NA} when the held-out cells contain a single
#'   class. Binary \code{actual} must be coded 0/1 and poisson \code{actual}
#'   must be non-negative integer counts, else an error is raised.
#' @examples
#' a <- matrix(rpois(100, 2), 10, 10); diag(a) <- NA
#' p <- a + 0.5
#' score_sir(a, p, "poisson")
#' @export
score_sir <- function(actual, predicted, family = c("poisson", "normal", "binomial"),
					  drop_diagonal = TRUE) {
	family <- match.arg(family)
	da <- dim(actual)
	dp <- dim(predicted)
	if (!is.null(da) || !is.null(dp)) {
		if (is.null(da) || is.null(dp) || !identical(da, dp)) {
			cli::cli_abort("{.arg actual} and {.arg predicted} must have identical dimensions.")
		}
	}
	a <- as.numeric(actual)
	p <- as.numeric(predicted)
	if (length(a) != length(p)) {
		cli::cli_abort("{.arg actual} and {.arg predicted} must have the same length.")
	}
	ok <- is.finite(a)

		# exclude one-mode self-ties for square arrays / stacks of square arrays
		if (isTRUE(drop_diagonal) && !is.null(da) && length(da) >= 2 && da[1] == da[2]) {
		if (length(da) == 2) {
			diagmask <- as.vector(!diag(TRUE, da[1]))
		} else {
			per <- prod(da[1:2])
			reps <- length(a) / per
			diagmask <- rep(as.vector(!diag(TRUE, da[1])), times = reps)
		}
		ok <- ok & diagmask
	}
	bad_pred <- ok & !is.finite(p)
	if (any(bad_pred)) {
		cli::cli_abort("{.arg predicted} has non-finite values on {sum(bad_pred)} scored cell{?s}.")
	}
	a <- a[ok]
	p <- p[ok]
	if (!length(a)) cli::cli_abort("No non-missing off-diagonal cells to score.")

		# validate the response domain so silently-wrong scores are not produced
			if (family == "poisson" && any(a < 0)) {
				cli::cli_abort("{.arg actual} has negative values; {.val poisson} scores require non-negative counts.")
			}
			if (family == "poisson" && any(abs(a - round(a)) > 1e-8)) {
				cli::cli_abort("{.arg actual} must be integer counts for {.val poisson} scoring.")
			}
			if (family == "poisson" && any(p < 0)) {
				cli::cli_abort("{.arg predicted} has negative values; {.val poisson} scores require non-negative means.")
			}
		if (family == "binomial" && !all(a %in% c(0, 1))) {
			cli::cli_abort("{.arg actual} must be 0/1 for {.val binomial} scoring.")
		}
		if (family == "binomial" && any(p < 0 | p > 1)) {
			cli::cli_abort("{.arg predicted} must be probabilities in [0, 1] for {.val binomial} scoring.")
		}

	rmse <- sqrt(mean((a - p)^2))
	mae  <- mean(abs(a - p))
	if (family == "normal") return(c(rmse = rmse, mae = mae))
	if (family == "poisson") {
		# unit deviance 2*(a*log(a/mu) - (a-mu)); the zero-count case is 2*mu.
		# Only positive observed counts need an epsilon guard when a predicted mean
		# is exactly zero. Using one mean per cell keeps the deviance non-negative.
		dev_terms <- numeric(length(a))
		pos <- a > 0
		mu_pos <- pmax(p[pos], .Machine$double.eps)
		dev_terms[pos] <- 2 * (a[pos] * log(a[pos] / mu_pos) - (a[pos] - mu_pos))
		dev_terms[!pos] <- 2 * p[!pos]
		dev <- mean(dev_terms)
		return(c(rmse = rmse, mae = mae, deviance = dev))
	}
	# binomial. logloss clamps predicted probabilities at eps, so hard 0/1
	# predictions (e.g. a last-value baseline) are charged ~-log(eps) per miss;
	# logloss/brier reward calibrated probabilities, auc ranks. prefer brier/auc
	# when comparing against hard-label baselines.
	eps <- 1e-15
	pc <- pmin(pmax(p, eps), 1 - eps)
	logloss <- -mean(a * log(pc) + (1 - a) * log(1 - pc))
	brier <- mean((a - p)^2)
	auc <- NA_real_
	# only a genuine two-class {0,1} target has a defined AUC
	if (setequal(unique(a), c(0, 1))) {
		rk <- rank(p)
		n_pos <- sum(a == 1)
		n_neg <- sum(a == 0)
		auc <- (sum(rk[a == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
	}
	c(rmse = rmse, mae = mae, logloss = logloss, brier = brier, auc = auc)
}

#' Rolling-origin cross-validation for a fitted sir model
#'
#' Performs expanding-window (rolling-origin) time-series cross-validation. For
#' each origin \code{o}, the model is refit on periods \code{1..o} and used to
#' forecast periods \code{o+1..o+horizon}, which are scored against the held-out
#' actuals. Training periods are always strictly before the test periods, but
#' forecasts are conditional on any future \code{Z} or dynamic-\code{W} slices
#' stored on \code{object}; interpret the result as conditional validation unless
#' those covariates are known, pre-specified, or generated as a scenario at the
#' forecast origin. A last-value-carried-forward naive baseline is scored
#' alongside the model so the result is a forecasting horse race.
#'
#' @param object a fitted \code{sir}/\code{sir_fit} object; its data, family, and
#'   structural settings are reused for the refits.
#' @param initial integer; length of the first training window. Default
#'   \code{ceiling(n_periods / 2)}.
#' @param horizon integer; forecast horizon scored at each origin (default 1, one
#'   step ahead).
#' @param origins optional integer vector of training-window end points; defaults
#'   to \code{seq(initial, n_periods - horizon)}. Pass a sparse vector (e.g.
#'   \code{seq(20, 90, by = 5)}) to subsample origins and cut cost.
#' @param baseline logical; also score the last-value-carried-forward baseline
#'   (default \code{TRUE}).
#' @param ... passed to the \code{\link{sir}} refit. The structural settings
#'   reused from \code{object} --- \code{family}, \code{method},
#'   \code{fix_receiver}, \code{symmetric}, \code{bipartite}, \code{calc_se},
#'   \code{W_recv} --- are reserved: passing any of them here is ignored with a
#'   warning (every fold must fit the same model as \code{object}).
#' @return a \code{sir_cv} object with the per-origin score table, the aggregate
#'   (mean) scores, the per-metric effective n (\code{eff_n}), the requested vs
#'   failed fold counts, and (if requested) the baseline aggregate.
#' @details
#' Each refit forces \code{calc_se = FALSE} (standard errors are not needed for
#' point forecasts and this is a large speed win). The cost is one model refit per
#' origin, and each refit's training window grows with the origin, so the default
#' (all origins from \code{initial} to \code{n_periods - horizon}) is roughly
#' quadratic in \code{n_periods} --- on the order of minutes for a hundred-period
#' series. For long series subsample \code{origins} (e.g.
#' \code{seq(initial, n_periods - horizon, by = 5)}). The model and the naive
#' baseline are scored on the same observed cells each origin so the comparison
#' is fair; non-finite model predictions on those cells make the fold fail rather
#' than silently changing the denominator. With \code{horizon > 1} the future
#' cells are pooled across the horizon into one score per origin. Folds that fail
#' to fit, forecast, or score are dropped with a warning and the aggregate
#' averages over those that succeeded. For a full-bilinear (\code{W_recv}) fit
#' the refits draw random restarts; the fit's \code{seed} is threaded through so
#' the result is reproducible.
#' @examples
#' dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
#' fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 1)
#' cv <- cv_sir(fit, initial = 12, origins = c(12, 15, 18))
#' cv
#' @export
cv_sir <- function(object, initial = NULL, horizon = 1L, origins = NULL,
				   baseline = TRUE, ...) {
	if (!inherits(object, "sir")) {
		cli::cli_abort(c(
			"{.arg object} must be a fitted {.cls sir} model, not {.obj_type_friendly {object}}.",
			"i" = "Fit first with {.fn sir}, then pass the result to {.fn cv_sir}."
		))
	}
	# horizon must be a positive integer (mirrors forecast()'s h guard)
	if (!is.numeric(horizon) || length(horizon) != 1 || !is.finite(horizon) ||
		horizon < 1 || horizon != round(horizon)) {
		cli::cli_abort("{.arg horizon} must be a single positive integer.")
	}
	horizon <- as.integer(horizon)

	Y <- object$Y
	W <- object$W
	X <- object$X
	Z <- object$Z
	family <- object$family
	Tt <- object$n_periods
	dyn <- isTRUE(object$dynamic_W)

	if (is.null(initial)) initial <- ceiling(Tt / 2)
	if (!is.numeric(initial) || length(initial) != 1L || !is.finite(initial) ||
		initial < 2 || initial != round(initial)) {
		cli::cli_abort("{.arg initial} must be a single integer >= 2.")
	}
	initial <- as.integer(initial)
	if (is.null(origins)) origins <- seq.int(initial, Tt - horizon)
	if (!is.numeric(origins) || !length(origins) || any(!is.finite(origins)) ||
		any(origins != round(origins))) {
		cli::cli_abort("{.arg origins} must be a non-empty integer vector.")
	}
	origins <- as.integer(origins)
	if (min(origins) < 2 || max(origins) + horizon > Tt) {
		cli::cli_abort(c(
			"Invalid {.arg origins}/{.arg initial}/{.arg horizon} for T = {Tt}.",
			"i" = "Each origin must satisfy 2 <= origin and origin + horizon <= {Tt}; the auto-default {.arg initial} may be too small for a short series."
		))
	}

	# structural settings reused on every refit (so CV fits the same model)
	struct <- list(
		fix_receiver = isTRUE(object$fix_receiver),
		symmetric    = isTRUE(object$symmetric),
		bipartite    = isTRUE(object$bipartite),
		method       = if (!is.null(object$method)) object$method else "ALS"
	)
	W_recv <- object$W_recv   # NULL unless a full-bilinear bipartite fit

	# the refit hard-codes these six; passing them via ... would collide in
	# do.call(sir). strip (with a warning) so a stray calc_se=/method= is a clear
	# message rather than an opaque per-fold failure.
	extra <- list(...)
	reserved <- c("family", "method", "fix_receiver", "symmetric", "bipartite",
				  "calc_se", "W_recv")
	clash <- intersect(names(extra), reserved)
	if (length(clash)) {
		cli::cli_warn(c(
			"Ignoring {.arg {clash}} passed via {.arg ...}: {.fn cv_sir} reuses the fitted model's settings on every refit.",
			"i" = "These are fixed to match {.arg object}: {.field {reserved}}."
		))
		extra <- extra[setdiff(names(extra), reserved)]
	}
	# thread the fit's seed so full-bilinear refits are reproducible, unless the
	# caller supplied one explicitly.
	if (is.null(extra$seed) && !is.null(object$seed)) extra$seed <- object$seed

	rows <- list()
	first_err <- NULL
	n_failed <- 0L
	for (o in origins) {
		Ytr <- Y[, , seq_len(o), drop = FALSE]
		Xtr <- X[, , seq_len(o), drop = FALSE]
		Ztr <- if (is.null(Z)) NULL else Z[, , , seq_len(o), drop = FALSE]
		Wtr <- if (dyn) W[, , , seq_len(o), drop = FALSE] else W

		refit_args <- c(list(Ytr, W = Wtr, X = Xtr, Z = Ztr, family = family,
						   method = struct$method, fix_receiver = struct$fix_receiver,
						   symmetric = struct$symmetric, bipartite = struct$bipartite,
						   calc_se = FALSE), extra)
		if (!is.null(W_recv)) refit_args$W_recv <- W_recv
		fit_o <- tryCatch(do.call(sir, refit_args),
						  error = function(e) { if (is.null(first_err)) first_err <<- conditionMessage(e); NULL })
		if (is.null(fit_o)) { n_failed <- n_failed + 1L; next }
		if (!isTRUE(fit_o$convergence)) {
			if (is.null(first_err)) first_err <- "A fold refit did not converge."
			n_failed <- n_failed + 1L
			next
		}

		Zf <- if (is.null(Z)) NULL else Z[, , , (o + 1):(o + horizon), drop = FALSE]
		Wf <- if (dyn) W[, , , (o + 1):(o + horizon), drop = FALSE] else NULL
		fc <- tryCatch(
			forecast(fit_o, h = horizon, Z_future = Zf, W_future = Wf),
			error = function(e) { if (is.null(first_err)) first_err <<- conditionMessage(e); NULL }
		)
		if (is.null(fc)) { n_failed <- n_failed + 1L; next }

		actual <- Y[, , (o + 1):(o + horizon), drop = FALSE]
		# Score the model and the naive baseline on the SAME observed cells. A
		# non-finite model prediction on those cells must fail the fold rather
		# than quietly shrinking the score denominator.
		mask <- is.finite(actual)
		base_pred <- NULL
		if (baseline) {
			base_pred <- array(Y[, , o], dim = dim(actual))
			mask <- mask & is.finite(base_pred)
		}
		actual_m <- actual; actual_m[!mask] <- NA
		fc_m <- fc; fc_m[!mask] <- NA
		drop_diagonal <- !isTRUE(object$bipartite)
		sc <- tryCatch(
			score_sir(actual_m, fc_m, family, drop_diagonal = drop_diagonal),
			error = function(e) { if (is.null(first_err)) first_err <<- conditionMessage(e); NULL }
		)
		if (is.null(sc)) { n_failed <- n_failed + 1L; next }
		row <- data.frame(origin = o, as.list(sc))
		if (baseline) {
			base_pred[!mask] <- NA
			bsc <- tryCatch(
				score_sir(actual_m, base_pred, family, drop_diagonal = drop_diagonal),
				error = function(e) { if (is.null(first_err)) first_err <<- conditionMessage(e); NULL }
			)
			if (is.null(bsc)) { n_failed <- n_failed + 1L; next }
			names(bsc) <- paste0("naive_", names(bsc))
			row <- cbind(row, as.list(bsc))
		}
		rows[[length(rows) + 1]] <- row
	}
	if (!length(rows)) {
		cli::cli_abort(c(
			"All CV folds failed to fit/forecast.",
			if (!is.null(first_err)) c("x" = "First error: {first_err}")
		))
	}
	if (n_failed > 0) {
		cli::cli_warn("{n_failed} of {length(origins)} CV fold{?s} failed and {?was/were} dropped; aggregates average over the {length(rows)} that succeeded.")
	}

	scores <- do.call(rbind, rows)
	score_cols <- setdiff(names(scores), "origin")
	aggregate <- colMeans(scores[, score_cols, drop = FALSE], na.rm = TRUE)
	# per-metric effective n (na.rm can drop origins for, e.g., single-class AUC)
	eff_n <- vapply(scores[, score_cols, drop = FALSE],
					function(col) sum(is.finite(col)), integer(1))
	structure(list(scores = scores, aggregate = aggregate, eff_n = eff_n,
				   family = family, n_origins = nrow(scores),
				   n_requested = length(origins), n_failed = n_failed,
				   horizon = horizon, initial = initial, call = match.call()),
			  class = "sir_cv")
}

#' Print rolling-origin cross-validation results
#'
#' @param x a \code{sir_cv} object returned by \code{\link{cv_sir}}.
#' @param digits integer; number of decimal places to display.
#' @param ... unused.
#' @return Invisibly returns \code{x}.
#' @method print sir_cv
#' @export
print.sir_cv <- function(x, digits = 3, ...) {
	cli::cli_h2("Rolling-origin cross-validation ({x$n_origins} origins, horizon {x$horizon})")
	cli::cli_text("Family: {.val {x$family}}")
	if (!is.null(x$n_failed) && x$n_failed > 0) {
		cli::cli_text("{.emph {x$n_failed} of {x$n_requested} requested fold{?s} failed and {?was/were} dropped.}")
	}
	ag <- x$aggregate
	model_cols <- ag[!grepl("^naive_", names(ag))]
	cli::cli_text("")
	cli::cli_text("{.strong Model (out-of-sample, averaged over origins):}")
	for (nm in names(model_cols)) {
		cli::cli_text("  {nm}: {.val {round(model_cols[[nm]], digits)}}")
	}
	naive_cols <- ag[grepl("^naive_", names(ag))]
	if (length(naive_cols)) {
		cli::cli_text("")
		cli::cli_text("{.strong Naive (last value carried forward):}")
		for (nm in names(naive_cols)) {
			cli::cli_text("  {sub('naive_', '', nm)}: {.val {round(naive_cols[[nm]], digits)}}")
		}
	}
	invisible(x)
}
