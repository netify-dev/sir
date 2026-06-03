# forecasting for fitted sir models: iterated plug-in out-of-sample prediction.
# reuses predict()/eta_tab for the bilinear algebra; the only new logic is the
# temporal recursion that turns a predicted mean into the next-period lag X.

#' @importFrom generics forecast
#' @importFrom stats predict
#' @export
generics::forecast

#' Forecast future networks from a fitted sir model
#'
#' Produces out-of-sample expected networks \code{h} steps beyond the end of the
#' estimation sample. Forecasts are iterated (plug-in): each period's predicted
#' mean is fed forward as the next-period lagged signal \code{X}, transformed the
#' same way \code{\link{sim_sir}} builds \code{X} from a lagged outcome
#' (\code{log(Y + 1) / infl_scale} for \code{poisson}, raw \code{Y / infl_scale}
#' otherwise). This assumes \code{X} was built with the same \code{infl_scale}
#' used below; if it was not, supply \code{infl_scale} --- the function warns
#' when the stored \code{X} disagrees. Forecasts are conditional on any supplied
#' future covariates, so \code{Z_future} and \code{W_future} should be values
#' that are known, pre-specified, or generated as a scenario at the forecast
#' origin.
#'
#' @param object a fitted \code{sir}/\code{sir_fit} object.
#' @param h integer >= 1; forecast horizon (number of future periods).
#' @param Z_future exogenous covariates for the future periods. Required when the
#'   model has \code{q > 0} predictors. These are conditioning values, not
#'   forecast by \code{forecast()}: pass only covariates known or fixed at the
#'   forecast origin, or values from an explicit scenario. Either an
#'   \code{n1 x n2 x q x h} array, or (when \code{h == 1}) an
#'   \code{n1 x n2 x q} array. Ignored with a message when \code{q == 0}.
#' @param W_future future influence covariates; required only when
#'   \code{object$dynamic_W} is \code{TRUE} and the model has influence
#'   covariates. As with \code{Z_future}, pass known, fixed, or scenario values.
#'   An \code{n1 x n1 x p x h} array (or \code{n1 x n1 x p} when \code{h == 1}).
#'   Static-W models reuse \code{object$W}.
#' @param Y_last optional \code{n1 x n2} matrix giving the most recent observed
#'   outcome that seeds the first forecast lag. Defaults to the last period of
#'   \code{object$Y}. Missing (NA) cells are treated as a zero lag.
#' @param infl_scale optional positive scalar; the divisor applied to the lagged
#'   outcome when building the forecast \code{X}. Defaults to \code{max(m - 1, 1)}
#'   for a one-mode network, \code{max(n1 - 1, 1)} for a sender-side bipartite
#'   fit, and \code{sqrt((n1 - 1)(n2 - 1))} for a full-bilinear bipartite fit.
#'   Override it when \code{X} was constructed with a different scaling.
#' @param ... unused.
#'
#' @return an \code{n1 x n2 x h} array of expected outcomes on the response
#'   scale, with actor dimnames carried from \code{object$Y} when present and a
#'   third dimension labelled \code{h1, h2, ...}. Carries attribute
#'   \code{"family"}.
#'
#' @details
#' For \code{h == 1} the forecast equals \code{predict(object, newdata =
#' list(X = X_next, ...))} where \code{X_next} is the transform of \code{Y_last}
#' (off the diagonal for a square one-mode network, whose self-ties are NA in the
#' forecast but computed by \code{predict}). For \code{h > 1} the period-\code{s}
#' predicted mean becomes the lag for period \code{s + 1}. This is a plug-in
#' point forecast: it feeds the conditional mean forward and does not propagate
#' forecast uncertainty, so it produces no predictive interval; for a nonlinear
#' link (\code{poisson}/\code{binomial}) the multi-step path is a biased
#' approximation of the true conditional mean (Jensen's inequality), growing with
#' the horizon. The diagonal (self-ties) is not a meaningful forecast for a
#' square one-mode network and is returned as \code{NA}.
#'
#' @examples
#' dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
#' train_t <- 19
#' fit <- sir(dat$Y[, , 1:train_t], W = dat$W,
#'            X = dat$X[, , 1:train_t],
#'            Z = dat$Z[, , , 1:train_t],
#'            family = "poisson", fix_receiver = TRUE, seed = 1)
#' # one-step-ahead forecast, conditional on pre-specified future Z.
#' # In this simulated example we reuse the held-out Z slice to show the API.
#' Z_next <- dat$Z[, , , train_t + 1, drop = FALSE]
#' fc1 <- forecast(fit, h = 1, Z_future = Z_next)
#' dim(fc1)
#' score_sir(dat$Y[, , train_t + 1, drop = FALSE], fc1, "poisson")
#' @export
forecast.sir_fit <- function(object, h = 1L, Z_future = NULL, W_future = NULL,
							 Y_last = NULL, infl_scale = NULL, ...) {
	if (!is.numeric(h) || length(h) != 1 || !is.finite(h) || h < 1 || h != round(h)) {
		cli::cli_abort("{.arg h} must be a single positive integer.")
	}
	h <- as.integer(h)
	family <- object$family
	m  <- object$m
	n1 <- object$n1
	n2 <- object$n2
	# the influence-carrying lag is divided by the number of partners summed
	# over. for square one-mode and sender-only bipartite fits this is one-sided;
	# full-bilinear bipartite fits use a symmetric two-sided scale. an explicit
	# infl_scale overrides this.
	bipartite <- isTRUE(object$bipartite)
	if (is.null(infl_scale)) {
		infl_scale <- if (bipartite && isTRUE(object$full_bilinear)) {
			sqrt(max(n1 - 1, 1) * max(n2 - 1, 1))
		} else if (bipartite) {
			max(n1 - 1, 1)
		} else {
			max(m - 1, 1)
		}
	} else if (!is.numeric(infl_scale) || length(infl_scale) != 1 || !is.finite(infl_scale) || infl_scale <= 0) {
		cli::cli_abort("{.arg infl_scale} must be a single positive number.")
	}
	q <- if (is.null(object$Z)) 0L else dim(object$Z)[3]
	dyn <- isTRUE(object$dynamic_W) && object$p > 0L

	# turn an outcome/mean into the next-period lagged signal X. this assumes the
	# same lag construction sim_sir uses; if the fit's stored X disagrees the
	# forecast would be mis-scaled, so warn when it does (see check below).
	to_lag <- function(Yprev) {
		x <- if (family == "poisson") log(Yprev + 1) else Yprev
		x[is.na(x)] <- 0
		x / infl_scale
	}

	# guard against a silently mis-scaled forecast: if the X the model was fit on
	# does not match to_lag(Y_prev) on the training periods, the lag convention
	# differs and the forecast scale is wrong. warn (do not abort) so a user with
	# a deliberately different X can still proceed via infl_scale.
	if (object$p > 0L && object$n_periods >= 2 &&
		!is.null(object$X) && !is.null(object$Y)) {
		t_chk <- object$n_periods
		implied <- to_lag(object$Y[, , t_chk - 1])
		stored  <- object$X[, , t_chk]
		off <- if (n1 == n2 && !bipartite) (row(implied) != col(implied)) else array(TRUE, dim(implied))
		denom <- max(stats::sd(stored[off]), 1e-8)
		rel <- max(abs(implied[off] - stored[off]), na.rm = TRUE) / denom
		if (is.finite(rel) && rel > 0.05) {
			cli::cli_warn(c(
				"The forecast lag transform does not match the {.field X} this model was fit on.",
				"i" = "{.fn forecast} assumes {.code X = log(Y+1)/infl_scale} (poisson) or {.code Y/infl_scale}; supply {.arg infl_scale} or rebuild {.field X} consistently."
			))
		}
	}

	# validate / normalise future covariates
	if (q > 0) {
		if (is.null(Z_future)) {
			cli::cli_abort("{.arg Z_future} is required when the model has q = {q} exogenous covariate{?s}.")
		}
		Z_future <- .sir_normalise_future(Z_future, n1, n2, q, h, "Z_future")
	} else if (!is.null(Z_future)) {
		cli::cli_inform("Model has no exogenous covariates (q = 0); ignoring {.arg Z_future}.")
	}
	if (dyn) {
		p <- dim(object$W)[3]
		if (is.null(W_future)) {
			cli::cli_abort("{.arg W_future} is required for a dynamic-W model.")
		}
		W_future <- .sir_normalise_future(W_future, n1, n1, p, h, "W_future")
	} else if (!is.null(W_future)) {
		cli::cli_inform("Model has static W; ignoring {.arg W_future}.")
	}

	# seed the first lag from the last observed outcome
	if (is.null(Y_last)) {
		Y_last <- object$Y[, , object$n_periods]
	} else if (!is.matrix(Y_last) || nrow(Y_last) != n1 || ncol(Y_last) != n2) {
		cli::cli_abort("{.arg Y_last} must be a {n1} x {n2} matrix.")
	}

	out <- array(NA_real_, dim = c(n1, n2, h))
	Yprev <- Y_last
		for (s in seq_len(h)) {
			X_s <- array(to_lag(Yprev), dim = c(n1, n2, 1))
			W_s <- if (dyn) array(W_future[, , , s], dim = c(n1, n1, dim(W_future)[3])) else object$W
			Z_s <- if (q > 0) array(Z_future[, , , s], dim = c(n1, n2, q, 1)) else NULL
			nd <- list(X = X_s, Z = Z_s)
			if (object$p > 0L) nd$W <- W_s
			mu_s <- predict(object, newdata = nd, type = "response")
		mu_s <- mu_s[, , 1]
		# self-ties are not modelled for a square one-mode network, so the
		# diagonal carries no forecast (it is NA in W) and is returned as NA.
		# a genuinely bipartite-square fit (distinct row/col populations) has
		# real diagonal cells, so leave them in place.
		if (n1 == n2 && !bipartite) diag(mu_s) <- NA
		out[, , s] <- mu_s
		Yprev <- mu_s
	}

	dn <- dimnames(object$Y)
	dimnames(out) <- list(if (!is.null(dn)) dn[[1]] else NULL,
						  if (!is.null(dn)) dn[[2]] else NULL,
						  paste0("h", seq_len(h)))
	attr(out, "family") <- family
	out
}

#' @rdname forecast.sir_fit
#' @export
forecast.sir <- forecast.sir_fit

# coerce a future-covariate array to n1 x n2 x k x h (3D allowed when h == 1)
.sir_normalise_future <- function(arr, d1, d2, k, h, nm) {
	da <- dim(arr)
	if (length(da) == 3 && h == 1 && da[1] == d1 && da[2] == d2 && da[3] == k) {
		return(array(arr, dim = c(d1, d2, k, 1)))
	}
	if (length(da) == 4 && da[1] == d1 && da[2] == d2 && da[3] == k && da[4] == h) {
		return(arr)
	}
	# a common slip: with k == 1, Z[,,,T] drops the singleton covariate axis to a
	# bare matrix. point the user at drop = FALSE.
	if (k == 1 && h == 1 && length(da) == 2 && da[1] == d1 && da[2] == d2) {
		cli::cli_abort(c(
			"{.arg {nm}} lost its covariate dimension (it is a {d1} x {d2} matrix, not {d1} x {d2} x 1).",
			"i" = "Subsetting like {.code Z[, , , T]} drops the singleton axis; use {.code drop = FALSE} or wrap in {.code array(., dim = c({d1}, {d2}, 1))}."
		))
	}
	cli::cli_abort("{.arg {nm}} must be {d1} x {d2} x {k} x {h} (or {d1} x {d2} x {k} when h = 1).")
}
