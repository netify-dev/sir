# re-export the broom/generics verbs so users can call tidy()/glance()/augment()
# directly after attaching sir, without also attaching broom or generics.

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics augment
#' @export
generics::augment

#' Tidy a SIR Model into a Data Frame of Coefficients
#'
#' Returns a data frame with one row per estimated parameter, in the layout
#' expected by \pkg{broom} consumers such
#' as \pkg{modelsummary} and \pkg{gtsummary}. Each parameter is tagged by its
#' role (\code{component}: exogenous \code{theta}, sender \code{alpha}, receiver
#' \code{beta}, or shared \code{gamma} for a symmetric fit).
#'
#' @param x A fitted \code{sir} object from \code{\link{sir}}.
#' @param conf.int Logical; if \code{TRUE}, add \code{conf.low}/\code{conf.high}
#'   Wald interval columns. Default \code{FALSE}.
#' @param conf.level Confidence level for the interval. Default 0.95.
#' @param se.type Which standard errors to report: \code{"cluster"} (default;
#'   actor-clustered sandwich, each cell scored onto both endpoint actors, for
#'   directed, symmetric, and dynamic (4D) \code{W} fits), \code{"classical"}
#'   (inverse-Hessian), or \code{"robust"} (HC0 sandwich). Ignored when the fit
#'   carries
#'   \code{se_source == "jackknife"} (analytic SEs were unavailable, so the
#'   delete-one-actor jackknife standard errors are reported for every type).
#' @param ... Unused, for generic compatibility.
#'
#' @return A data frame with columns \code{term}, \code{component},
#'   \code{estimate}, \code{std.error}, \code{statistic}, \code{p.value}, and
#'   (optionally) \code{conf.low}, \code{conf.high}.
#'
#' @examples
#' dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
#' fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 1)
#' tidy(fit)
#' tidy(fit, conf.int = TRUE)
#' @importFrom generics tidy
#' @export
tidy.sir <- function(x, conf.int = FALSE, conf.level = 0.95,
					 se.type = c("cluster", "classical", "robust"), ...) {
	se.type <- match.arg(se.type)
	summ <- x$summ
	term <- rownames(summ)
	if (is.null(term)) term <- paste0("p", seq_len(nrow(summ)))

	# map the name tag to a component; ^\(betaW catches both (betaW) and (betaWr)
	component <- ifelse(grepl("^\\(alphaW", term), "alpha",
				 ifelse(grepl("^\\(betaW", term), "beta",
				 ifelse(grepl("^\\(gammaW", term), "gamma", "theta")))

	# se + reference distribution mirror confint.sir: cluster sandwich SEs use a
	# t(G-1) reference (df on the vcov "cluster_df" attr), classical/HC0 use normal
	is_sym <- isTRUE(x$symmetric) && identical(x$operator, "symmetric")
	crit_df <- Inf
	se <- if (identical(x$se_source, "jackknife")) {
		# analytic SEs were unavailable, so sir() attached the delete-one-actor
		# jackknife covariance; vcov()/confint() use it for every type, so does tidy
		sqrt(pmax(diag(x$vcov), 0))
	} else if (se.type == "classical") {
		if (isFALSE(x$se_reliable)) rep(NA_real_, nrow(summ)) else summ$se
	} else if (se.type == "robust" && !is_sym) {
		if (isFALSE(x$se_reliable)) rep(NA_real_, nrow(summ)) else summ$rse
	} else {
		Vc <- vcov(x, type = se.type)
		crit_df <- attr(Vc, "cluster_df")
		sqrt(diag(Vc))
	}
	tstat <- if (!is.null(se)) summ$coef / se else rep(NA_real_, nrow(summ))
	pval <- if (is.null(se)) {
		rep(NA_real_, nrow(summ))
	} else if (is.null(crit_df) || !is.finite(crit_df) || crit_df < 1) {
		2 * stats::pnorm(abs(tstat), lower.tail = FALSE)
	} else {
		2 * stats::pt(abs(tstat), df = crit_df, lower.tail = FALSE)
	}

	# flag the silent all-NA-inference case so it doesn't land unnoticed in a
	# publication table (happens when the fit used calc_se = FALSE)
	if (is.null(se) || all(is.na(se))) {
		cli::cli_inform(c(
			"i" = "Standard errors are unavailable, so {.field std.error}/{.field statistic}/{.field p.value} are NA.",
			" " = "Refit with {.code calc_se = TRUE}, use an explicit supported {.arg se.type}, or use {.fn boot_sir} for inference."
		))
	} else if (isFALSE(x$se_reliable) && se.type %in% c("classical", "robust")) {
		cli::cli_inform(c(
			"i" = "Hessian-based SEs were marked unreliable; classical/robust inference columns are NA.",
			" " = "Use {.code boot_sir()} or refit/simplify the model before reporting inference."
		))
	}

	out <- data.frame(
		term = term,
		component = component,
		estimate = summ$coef,
		std.error = if (is.null(se)) NA_real_ else se,
		statistic = tstat,
		p.value = pval,
		stringsAsFactors = FALSE,
		row.names = NULL
	)

	if (conf.int) {
		ci <- confint(x, level = conf.level, se.type = se.type)
		out$conf.low  <- ci[, 1]
		out$conf.high <- ci[, 2]
	}

	out
}

#' One-Row Model Summary for a SIR Model
#'
#' Returns a single-row data frame of model-level statistics, in the layout
#' \pkg{broom} consumers expect for model comparison tables.
#'
#' @param x A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Unused, for generic compatibility.
#'
#' @return A one-row data frame with columns \code{nobs}, \code{df}
#'   (number of estimated parameters), \code{logLik}, \code{AIC}, \code{BIC},
#'   \code{n_nodes}, \code{n_periods}, \code{n_influence_covar} (p),
#'   \code{family}, \code{method}, and \code{converged}.
#'
#' @examples
#' dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
#' fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 1)
#' glance(fit)
#' @importFrom generics glance
#' @export
glance.sir <- function(x, ...) {
	ll <- tryCatch(as.numeric(logLik(x)), error = function(e) NA_real_)
	df <- tryCatch(attr(logLik(x), "df"), error = function(e) length(x$tab))
	data.frame(
		nobs = x$nobs,
		df = df,
		logLik = ll,
		AIC = tryCatch(AIC(x), error = function(e) NA_real_),
		BIC = tryCatch(BIC(x), error = function(e) NA_real_),
		# network-size metadata, natural for a relational model
		n_nodes = if (!is.null(x$m)) x$m else x$n1,
		n_periods = x$n_periods,
		# sender-side p, plus receiver-side p2 for a full-bilinear bipartite fit
		n_influence_covar = x$p + (if (!is.null(x$p2)) x$p2 else 0L),
		family = x$family,
		method = x$method,
		converged = isTRUE(x$convergence),
		stringsAsFactors = FALSE,
		row.names = NULL
	)
}

#' Augment Data With SIR Model Fitted Values and Residuals
#'
#' Returns the long-format dyad-time table of observed outcomes, fitted values,
#' and residuals from a fitted SIR model. Because SIR data are arrays rather
#' than a single data frame, the returned table is built from the model's own
#' \code{Y}/\code{fitted.values}/\code{residuals} with sender, receiver, and
#' time index columns.
#'
#' @param x A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Unused, for generic compatibility.
#'
#' @return A data frame with columns \code{sender}, \code{receiver},
#'   \code{time}, \code{.observed}, \code{.fitted}, \code{.resid} (response
#'   residual), and \code{.resid_pearson}. Diagonal (self-tie) cells and any
#'   cells excluded from the likelihood are dropped.
#'
#' @examples
#' dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
#' fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 1)
#' head(augment(fit))
#' @importFrom generics augment
#' @export
augment.sir <- function(x, ...) {
	Y <- x$Y
	fit <- x$fitted.values
	dn <- dimnames(Y)
	d <- dim(Y)
	idx <- expand.grid(sender = seq_len(d[1]), receiver = seq_len(d[2]),
					   time = seq_len(d[3]))
	# use labels when available
	lab <- function(v, k) if (!is.null(dn) && !is.null(dn[[k]])) dn[[k]][v] else v
	out <- data.frame(
		sender = lab(idx$sender, 1),
		receiver = lab(idx$receiver, 2),
		time = lab(idx$time, 3),
		.observed = as.vector(Y),
		.fitted = as.vector(fit),
		.resid = if (!is.null(x$residuals$response)) as.vector(x$residuals$response) else as.vector(Y) - as.vector(fit),
		.resid_pearson = if (!is.null(x$residuals$pearson)) as.vector(x$residuals$pearson) else NA_real_,
		stringsAsFactors = FALSE
	)
	# drop cells not in the likelihood (NA observed, e.g. diagonal / masked)
	out[!is.na(out$.observed), , drop = FALSE]
}
