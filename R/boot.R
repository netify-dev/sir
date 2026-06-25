
# delete-one-actor jackknife covariance for the influence parameters. deleting
# (rather than resampling-with-replacement) avoids duplicating actors into the
# bilinear A X B' sums, which would attenuate the influence coefficients. for a
# square directed network the same actor is dropped from both axes together; for
# a bipartite network senders and receivers are dropped separately and the two
# one-way jackknife covariances are summed as a delete-one-margin sensitivity
# check.
.sir_dyad_jackknife <- function(sir_fit, trace = FALSE) {
	Y <- sir_fit$Y; W <- sir_fit$W; X <- sir_fit$X; Z <- sir_fit$Z
	W_recv <- sir_fit$W_recv
	family <- sir_fit$family
	n1 <- dim(Y)[1]; n2 <- dim(Y)[2]
	P <- length(sir_fit$tab)
	bipartite <- isTRUE(sir_fit$bipartite)
	fr <- isTRUE(sir_fit$fix_receiver)
	sym <- isTRUE(sir_fit$symmetric)
	dynamic_W <- isTRUE(sir_fit$dynamic_W)

	z_sub <- function(i1, i2) {
		if (is.null(Z)) return(NULL)
		if (length(dim(Z)) == 4) Z[i1, i2, , , drop = FALSE] else Z[i1, i2, , drop = FALSE]
	}
	# refit on the sub-network induced by keeping sender set i1, receiver set i2
	refit <- function(i1, i2) {
		Y_s <- Y[i1, i2, , drop = FALSE]
		if (!bipartite && length(i1) == length(i2) && all(i1 == i2)) {
			for (tt in seq_len(dim(Y_s)[3])) diag(Y_s[, , tt]) <- NA
		}
		args <- list(Y_s,
					 W = if (dynamic_W) W[i1, i1, , , drop = FALSE] else W[i1, i1, , drop = FALSE],
					 X = X[i1, i2, , drop = FALSE], Z = z_sub(i1, i2),
					 family = family, method = "ALS", calc_se = FALSE,
					 bipartite = bipartite, symmetric = sym,
					 fix_receiver = fr, max_iter = 100, tol = 1e-6)
		if (!is.null(W_recv)) args$W_recv <- W_recv[i2, i2, , drop = FALSE]
		# per-refit warnings (small-subnetwork gain, etc.) are expected across
		# delete-one fits and would flood the console, so silence them here
		tab <- tryCatch(suppressWarnings({
			fit <- do.call(sir, args)
			if (!isTRUE(fit$convergence)) return(NULL)
			fit$tab
		}), error = function(e) NULL)
		if (is.null(tab) || length(tab) != P) NULL else tab
	}
	# symmetric A is identified up to global sign; flip each refit's gamma block to
	# agree with the point estimate so resampled signs do not cancel in the variance
	align_sign <- function(mat) {
		if (!sym || is.null(mat) || !nrow(mat)) return(mat)
		q <- if (is.null(sir_fit$q)) 0L else sir_fit$q
		p <- if (is.null(sir_fit$p)) 0L else sir_fit$p
		gidx <- q + seq_len(p)
		if (ncol(mat) < max(gidx)) return(mat)
		pg <- sir_fit$tab[gidx]
		flip <- vapply(seq_len(nrow(mat)),
					   function(k) sum(mat[k, gidx] * pg, na.rm = TRUE) < 0, logical(1))
		mat[flip, gidx] <- -mat[flip, gidx]
		mat
	}

	# jackknife covariance from the delete-one estimates of one margin
	jk_cov <- function(mat) {
		g <- nrow(mat)
		if (is.null(g) || g < 2) return(matrix(0, P, P))
		mbar <- colMeans(mat)
		(g - 1) / g * Reduce(`+`, lapply(seq_len(g), function(k) {
			d <- mat[k, ] - mbar; outer(d, d)
		}))
	}

		if (bipartite || n1 != n2) {
			# bipartite: drop senders and receivers separately, sum the covariances
			Ms <- align_sign(do.call(rbind, Filter(Negate(is.null),
				lapply(seq_len(n1), function(a) refit(setdiff(seq_len(n1), a), seq_len(n2))))))
			Mr <- align_sign(do.call(rbind, Filter(Negate(is.null),
				lapply(seq_len(n2), function(a) refit(seq_len(n1), setdiff(seq_len(n2), a))))))
			gs <- if (is.null(Ms)) 0L else nrow(Ms)
			gr <- if (is.null(Mr)) 0L else nrow(Mr)
			n_valid <- gs + gr
			n_total <- n1 + n2
			V <- jk_cov(Ms) + jk_cov(Mr)
			jack_est <- rbind(Ms, Mr)
			if (gs < 2L || gr < 2L) {
				cli::cli_warn("Dyad jackknife has fewer than two valid refits in at least one bipartite margin; intervals are a weak sensitivity check.")
			}
		} else {
			# square directed: drop the same actor from both axes together
			M <- align_sign(do.call(rbind, Filter(Negate(is.null),
				lapply(seq_len(n1), function(a) { k <- setdiff(seq_len(n1), a); refit(k, k) }))))
			n_valid <- if (is.null(M)) 0L else nrow(M)
			n_total <- n1
			V <- jk_cov(M)
			jack_est <- M
		}
		if (n_valid < 2) cli::cli_abort("Dyad jackknife produced too few valid refits.")
		if (n_valid < n_total) {
			cli::cli_warn("{n_total - n_valid} of {n_total} dyad jackknife refit{?s} failed or did not converge and {?was/were} dropped.")
		}

		se <- sqrt(pmax(diag(V), 0))
		point <- sir_fit$tab
		z <- stats::qnorm(0.975)
		list(se = se, ci_lo = point - z * se, ci_hi = point + z * se,
			 cov = V, jack_est = jack_est, n_valid = n_valid, n_total = n_total)
}

#' Bootstrap Inference for SIR Model Parameters
#'
#' Computes bootstrap standard errors and confidence intervals for SIR model
#' parameters. This is the recommended approach for inference when the Hessian
#' is singular or ill-conditioned, which is common in models with bilinear
#' influence terms (i.e., when \code{fix_receiver = FALSE}).
#'
#' Three bootstrap strategies are available:
#' \describe{
#'   \item{block}{Resamples whole time periods with replacement as independent
#'     period blocks. This preserves the within-period network dependence
#'     structure but not serial dependence across adjacent periods. Best when T
#'     is moderately large (T >= 10).}
#'   \item{dyad}{A delete-one-actor jackknife: each actor is dropped in turn and
#'     the model is refit on the induced sub-network, capturing the
#'     cross-sectional dyadic dependence (shared-sender / shared-receiver
#'     correlation) the block bootstrap ignores. Deletion (rather than
#'     resampling with replacement) avoids duplicating actors into the bilinear
#'     \eqn{A X B'} sums, which would attenuate the influence coefficients toward
#'     zero. For a square directed network the same actor is dropped from both
#'     axes together; for a bipartite network senders and receivers are dropped
#'     separately and the two one-way jackknife covariances are summed. Standard
#'     errors come from the jackknife covariance and intervals are normal
#'     (\code{estimate +/- z * se}). This is also the estimator \code{\link{sir}}
#'     reuses automatically when \code{calc_se = TRUE} cannot form analytic SEs (a
#'     singular/ill-conditioned Hessian or a full-bilinear bipartite fit): it
#'     attaches the same \code{$cov} so \code{vcov}/\code{confint}/\code{tidy}
#'     return jackknife inference without an explicit \code{boot_sir} call.}
#'   \item{parametric}{Simulates new outcome arrays from the fitted model
#'     using the estimated parameters and the specified family distribution.
#'     Better when T is small but the model is well-specified.}
#' }
#'
#' For \code{block} and \code{parametric}, each replicate refits the full SIR
#' model; replicates that fail to converge are dropped and reported, standard
#' errors are the column standard deviations of the successful replicates, and
#' confidence intervals use the percentile method. The \code{dyad} jackknife
#' instead enumerates all delete-one-actor refits and reports the jackknife
#' standard error and a normal interval.
#'
#' @param sir_fit A fitted \code{sir} object from \code{\link{sir}}.
#' @param R Integer. Number of bootstrap replicates. Default is 200. Increase
#'   to 500-1000 for publication-quality intervals.
#' @param type Character. Inference type: \code{"block"} (default) resamples
#'   whole time periods as independent blocks, preserving within-period network
#'   dependence but not serial dependence; \code{"dyad"} is a delete-one-actor
#'   jackknife on the induced sub-network (captures dyadic dependence, and the
#'   only inference path for full-bilinear bipartite fits); \code{"parametric"}
#'   simulates new outcomes from the fitted model. \code{R} is ignored for
#'   \code{"dyad"} (it enumerates all actor deletions).
#' @param seed Optional integer for reproducibility. Sets the random seed
#'   before resampling. For exact reproducibility when \code{cores > 1}, the
#'   "L'Ecuyer-CMRG" RNG is used so parallel workers draw independent streams;
#'   results are then reproducible given the same \code{seed} and \code{cores}.
#'   Serial runs (\code{cores = 1}) with the same \code{seed} are always
#'   reproducible.
#' @param trace Logical. If TRUE, shows progress and prints a verbose message
#'   every 10 serial replicates. Ignored when \code{cores > 1}.
#' @param cores Integer. Number of CPU cores to use. Default is 1 (serial). When
#'   greater than 1, replicates are run in parallel with
#'   \code{parallel::mclapply} (forking). Forking is unavailable on Windows, so
#'   there \code{cores > 1} falls back to serial with a one-time warning.
#'
#' @return An object of class \code{"boot_sir"} with components:
#' \describe{
#'   \item{coefs}{R x n_params matrix of bootstrap coefficient estimates.
#'     Rows for failed replicates contain NA.}
#'   \item{se}{Named numeric vector of bootstrap standard errors (one per
#'     parameter).}
#'   \item{cov}{For \code{type = "dyad"}, the jackknife variance-covariance
#'     matrix (NULL for the block/parametric bootstraps).}
#'   \item{ci_lo}{Lower 2.5\% percentile bounds.}
#'   \item{ci_hi}{Upper 97.5\% percentile bounds.}
#'   \item{point_est}{Point estimates from the original fit.}
#'   \item{param_names}{Character vector of parameter names.}
#'   \item{n_valid}{Number of successful bootstrap replicates.}
#'   \item{n_total}{Total number of replicates attempted.}
#'   \item{type}{The bootstrap type used.}
#'   \item{family}{The distribution family.}
#' }
#'
#' @examples
#' \donttest{
#' dat <- sim_sir(m = 8, T_len = 12, p = 2, q = 1, family = "poisson", seed = 1)
#' model <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z,
#'              family = "poisson", seed = 1)
#'
#' # use a larger R (for example 500+) for publication-quality intervals
#' boot_result <- boot_sir(model, R = 10, seed = 42)
#' print(boot_result)
#'
#' # use bootstrap CIs with confint
#' confint(model, boot = boot_result)
#' }
#'
#' @seealso \code{\link{confint.sir}} to use bootstrap intervals,
#'   \code{\link{confint.boot_sir}} for direct interval extraction.
#' @export
boot_sir <- function(sir_fit, R = 200, type = c("block", "parametric", "dyad"),
					 seed = NULL, trace = FALSE, cores = 1L) {
	type <- match.arg(type)
	if (!is.numeric(R) || length(R) != 1L || !is.finite(R) || R < 2 || R != round(R)) {
		cli::cli_abort("{.arg R} must be a single integer >= 2.")
	}
	R <- as.integer(R)
	if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
		cores < 1 || cores != round(cores)) {
		cli::cli_abort("{.arg cores} must be a positive integer.")
	}
	cores <- as.integer(cores)

	# inference on a fit that itself did not converge is unreliable; surface it.
	if (!is.null(sir_fit$convergence) && !isTRUE(sir_fit$convergence)) {
		cli::cli_warn(c(
			"The fitted model did not converge; bootstrap/jackknife inference on it may be unreliable.",
			"i" = "Refit until convergence (e.g. raise {.arg max_iter} or {.arg n_restarts}) before trusting these intervals."
		))
	}

	# parallel forking is unavailable on windows; fall back to serial
	if (cores > 1 && .Platform$OS.type == "windows") {
		cli::cli_warn(c(
			"Parallel bootstrap is unavailable on Windows (no forking).",
			"i" = "Falling back to serial ({.code cores = 1})."
		))
		cores <- 1L
	}

	# for reproducible parallel streams, switch to L'Ecuyer-CMRG before seeding.
	# restore the user's RNGkind and global RNG stream on exit so a seeded run
	# leaves no side effect (mirrors sir()/sim_sir()).
	if (!is.null(seed)) {
		if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
			old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
			on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
		} else {
			on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
				rm(".Random.seed", envir = globalenv()), add = TRUE)
		}
		if (cores > 1) {
			old_kind <- RNGkind("L'Ecuyer-CMRG")
			on.exit(RNGkind(old_kind[1]), add = TRUE)
		}
		set.seed(seed)
	}

	Y <- sir_fit$Y
	W <- sir_fit$W
	X <- sir_fit$X
	Z <- sir_fit$Z
	W_recv <- sir_fit$W_recv   # non-NULL only for full-bilinear bipartite fits
	family <- sir_fit$family
	T_len <- dim(Y)[3]
	n1 <- dim(Y)[1]
	n2 <- dim(Y)[2]
	n_params <- length(sir_fit$tab)
	p <- if (is.null(W)) 0L else dim(W)[3]
	q <- if (is.null(Z)) 0L else dim(Z)[3]

	dynamic_W <- isTRUE(sir_fit$dynamic_W)

	pnames <- rownames(sir_fit$summ)
	point_est <- sir_fit$tab
	names(point_est) <- pnames   # keep point_est names consistent with se / ci

	# dyad inference uses a delete-one-actor jackknife, not a with-replacement
	# resample: duplicating an actor inserts spurious zero self-pairs into the
	# bilinear A X B' sums and attenuates the influence coefficients toward zero.
	if (type == "dyad") {
		jk <- .sir_dyad_jackknife(sir_fit, trace = trace)
		coefs <- jk$jack_est
		colnames(coefs) <- pnames
		se <- jk$se; names(se) <- pnames
		ci_lo <- jk$ci_lo; ci_hi <- jk$ci_hi
		names(ci_lo) <- names(ci_hi) <- pnames
		cov <- jk$cov; dimnames(cov) <- list(pnames, pnames)
			result <- list(
				coefs = coefs, se = se, ci_lo = ci_lo, ci_hi = ci_hi, cov = cov,
				point_est = point_est, param_names = pnames,
				n_valid = jk$n_valid, n_total = jk$n_total,
				type = "dyad", family = family, interval = "normal-jackknife"
			)
		class(result) <- "boot_sir"
		return(result)
	}

	# one bootstrap replicate: resample/simulate and refit.
	# returns the public normalized coefficient vector, or NA on failure.
	run_rep <- function(b) {
		W_recv_b <- W_recv
		if (type == "block") {
			# resample time periods with replacement
			t_idx <- sample(1:T_len, T_len, replace = TRUE)
			Y_b <- Y[,, t_idx, drop = FALSE]
			X_b <- X[,, t_idx, drop = FALSE]
			Z_b <- if (!is.null(Z) && length(dim(Z)) == 4) {
				Z[,,, t_idx, drop = FALSE]
			} else {
				Z
			}
			# resample W's time dimension for dynamic (4D) W
			W_b <- if (dynamic_W) W[,,, t_idx, drop = FALSE] else W
		} else {
			# parametric: simulate from fitted model. symmetric fits use the
			# quadratic-form predictor; full-bilinear bipartite fits need the
			# two-sided linear predictor (B is not the identity).
			eta <- if (isTRUE(sir_fit$symmetric)) {
				eta_tab_symmetric(sir_fit$tab, W, X, Z, sir_fit$p, sir_fit$q)
			} else if (!is.null(W_recv)) {
				eta_tab_bipartite(sir_fit$tab, W, W_recv, X, Z,
								  sir_fit$p, sir_fit$p2, sir_fit$q)
			} else {
				eta_tab(sir_fit$tab, W, X, Z,
						fix_receiver = isTRUE(sir_fit$fix_receiver))
			}
			if (family == "normal") {
				sigma <- sqrt(sir_fit$sigma2)
				Y_b <- array(rnorm(length(eta), mean = eta, sd = sigma),
							 dim = dim(Y))
			} else if (family == "poisson") {
				lambda <- exp(eta)
				lambda[lambda > 1e6] <- 1e6
				Y_b <- array(rpois(length(eta), lambda = lambda),
							 dim = dim(Y))
			} else {
				prob <- 1 / (1 + exp(-eta))
				Y_b <- array(rbinom(length(eta), 1, prob),
							 dim = dim(Y))
			}
			# Preserve the exact estimation mask. This keeps off-diagonal missingness
			# missing and avoids dropping diagonal cells in square bipartite fits.
			Y_b[is.na(sir_fit$Y)] <- NA
			X_b <- X
			Z_b <- Z
			W_b <- W
		}

			tryCatch({
				refit_args <- list(Y_b, W_b, X_b, Z_b, family = family,
								   method = "ALS", calc_se = FALSE,
								   symmetric = isTRUE(sir_fit$symmetric),
								   fix_receiver = isTRUE(sir_fit$fix_receiver),
								   bipartite = isTRUE(sir_fit$bipartite),
								   kron_mode = isTRUE(sir_fit$kron_mode),
								   max_iter = 100, tol = 1e-6)
			# route to the full-bilinear estimator when the original fit used it
			if (!is.null(W_recv)) refit_args$W_recv <- W_recv_b
			fit_b <- do.call(sir, refit_args)
			if (!isTRUE(fit_b$convergence)) return(rep(NA_real_, n_params))
					fit_b$tab
			}, error = function(e) {
			# failed replicate
			rep(NA_real_, n_params)
		})
	}

		if (cores > 1) {
			# parallel path: mc.set.seed=TRUE makes each fork advance the
			# L'Ecuyer-CMRG stream, so streams are independent and reproducible.
			if (trace) cli::cli_inform("Bootstrapping {.val {R}} replicates on {.val {cores}} cores ...")
			rep_list <- parallel::mclapply(seq_len(R), run_rep,
										   mc.cores = cores, mc.set.seed = TRUE)
		} else {
			# serial path: trace controls both the progress bar and periodic messages
			rep_list <- vector("list", R)
			if (trace) cli::cli_progress_bar("Bootstrapping", total = R, clear = FALSE)
			for (b in seq_len(R)) {
				if (trace && b %% 10 == 0) cli::cli_inform("Bootstrap {.val {b}}/{.val {R}}")
				rep_list[[b]] <- run_rep(b)
				if (trace) cli::cli_progress_update()
			}
			if (trace) cli::cli_progress_done()
		}

	# assemble replicate coefficients into the R x n_params matrix.
	# guard against worker errors that yield non-conformable results.
	boot_coefs <- matrix(NA, R, n_params)
	colnames(boot_coefs) <- rownames(sir_fit$summ)
	for (b in seq_len(R)) {
		v <- rep_list[[b]]
		if (is.numeric(v) && length(v) == n_params) boot_coefs[b, ] <- v
	}

	# compute statistics from successful replicates
	valid <- apply(boot_coefs, 1, function(r) !any(is.na(r)))
	n_valid <- sum(valid)

	if (n_valid < 10) {
		cli::cli_warn("Only {.val {n_valid}} valid bootstrap replicates (of {.val {R}}). Results unreliable.")
	}
	if (n_valid < 2) {
		cli::cli_abort("Bootstrap produced fewer than 2 valid replicates; intervals and SEs are undefined.")
	}

	# symmetric A is identified up to global sign; align each replicate's gamma
	# block to the point estimate so resampled signs do not cancel
	if (isTRUE(sir_fit$symmetric)) {
		q <- if (is.null(sir_fit$q)) 0L else sir_fit$q
		p <- if (is.null(sir_fit$p)) 0L else sir_fit$p
		gidx <- q + seq_len(p)
		if (ncol(boot_coefs) >= max(gidx)) {
			pg <- point_est[gidx]
			for (b in which(valid)) {
				if (sum(boot_coefs[b, gidx] * pg, na.rm = TRUE) < 0) {
					boot_coefs[b, gidx] <- -boot_coefs[b, gidx]
				}
			}
		}
	}

	boot_se <- apply(boot_coefs[valid, , drop = FALSE], 2, sd)
	boot_ci <- apply(boot_coefs[valid, , drop = FALSE], 2,
					 quantile, probs = c(0.025, 0.975))

	result <- list(
		coefs = boot_coefs,
		se = boot_se,
		ci_lo = boot_ci[1, ],
		ci_hi = boot_ci[2, ],
		point_est = point_est,
		param_names = pnames,
		n_valid = n_valid,
		n_total = R,
		type = type,
		family = sir_fit$family,
		interval = "percentile"
	)
	class(result) <- "boot_sir"
	result
}

#' Print Bootstrap SIR Results
#'
#' Displays a table of point estimates, bootstrap standard errors, and 95\%
#' percentile confidence intervals.
#'
#' @param x A \code{boot_sir} object from \code{\link{boot_sir}}.
#' @param digits Number of significant digits. Default uses
#'   \code{getOption("digits") - 3}.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the \code{boot_sir} object.
#' @export
print.boot_sir <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
	cli::cli_text("\n")
	cli::cli_text("{.strong Bootstrap SIR Results}")
	cli::cli_text("Type: {.field {x$type}} | Replicates: {.val {x$n_valid}}/{.val {x$n_total}} valid")

	# build results table
	se_label <- if (identical(x$type, "dyad")) "Jackknife SE" else "Boot SE"
	tab <- cbind(
		Estimate = x$point_est,
		SE = x$se,
		`2.5 %` = x$ci_lo,
		`97.5 %` = x$ci_hi
	)
	colnames(tab)[2] <- se_label
	rownames(tab) <- x$param_names

	cli::cli_text("")
	print(round(tab, digits))
	cli::cli_text("")
	invisible(x)
}

#' Summary of Bootstrap SIR Results
#'
#' Displays detailed bootstrap results including the coefficient table,
#' significance indicators (whether the 95\% CI excludes zero), and
#' bootstrap distribution summaries (mean, sd, median).
#'
#' @param object A \code{boot_sir} object from \code{\link{boot_sir}}.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the \code{boot_sir} object.
#' @export
summary.boot_sir <- function(object, ...) {
	cli::cli_h1("Bootstrap SIR Results")
	cli::cli_text("Type: {.field {object$type}} | Family: {.field {object$family}}")
	cli::cli_text("Replicates: {.val {object$n_valid}} valid of {.val {object$n_total}} total ({.val {round(100 * object$n_valid / object$n_total, 1)}}%)")

	cli::cli_rule()
	cli::cli_h2("Parameter Estimates")

	tab <- cbind(
		Estimate = object$point_est,
		SE = object$se,
		`2.5 %` = object$ci_lo,
		`97.5 %` = object$ci_hi
	)
	colnames(tab)[2] <- if (identical(object$type, "dyad")) "Jackknife SE" else "Boot SE"
	rownames(tab) <- object$param_names

	# flag parameters where CI includes zero
	covers_zero <- object$ci_lo <= 0 & object$ci_hi >= 0
	sig <- ifelse(covers_zero, " ", "*")
	tab_print <- cbind(format(round(tab, 4), width = 10), ` ` = sig)

	print(noquote(tab_print))
	cli::cli_text("{.emph * = 95% CI excludes zero}")

	# distribution summary. for block/parametric this is the replicate
	# distribution; the dyad jackknife yields delete-one estimates, not a
	# sampling distribution, so its se comes from the jackknife covariance.
	valid <- apply(object$coefs, 1, function(r) !any(is.na(r)))
	boot_valid <- object$coefs[valid, , drop = FALSE]

	cli::cli_rule()
	if (identical(object$type, "dyad")) {
		cli::cli_h2("Delete-one-actor jackknife estimates")
	} else {
		cli::cli_h2("Bootstrap Distribution")
	}
	dist_tab <- data.frame(
		mean = colMeans(boot_valid),
		sd = apply(boot_valid, 2, sd),
		median = apply(boot_valid, 2, median),
		row.names = object$param_names
	)
	print(round(dist_tab, 4))

	invisible(object)
}

#' Confidence Intervals from Bootstrap SIR Results
#'
#' Computes confidence intervals from a \code{boot_sir} object. For \code{block}
#' and \code{parametric} bootstraps these are percentile intervals from the
#' replicate distribution. For the \code{dyad} jackknife they are normal
#' (\code{estimate +/- z * se}) intervals from the jackknife standard error;
#' the jackknife produces a covariance, not a replicate distribution, so
#' percentiles do not apply.
#'
#' @param object A \code{boot_sir} object from \code{\link{boot_sir}}.
#' @param parm Character vector of parameter names or integer indices. If
#'   NULL (default), returns intervals for all parameters.
#' @param level Confidence level between 0 and 1. Default is 0.95.
#' @param ... Additional arguments (unused).
#' @return A matrix with one row per parameter and columns for the lower
#'   and upper bounds, labeled by percentage.
#' @export
confint.boot_sir <- function(object, parm = NULL, level = 0.95, ...) {
	if (!inherits(object, "boot_sir")) {
		cli::cli_abort("{.arg object} must be a {.cls boot_sir} object.")
	}
	if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
		level <= 0 || level >= 1) {
		cli::cli_abort("{.arg level} must be a single number between 0 and 1.")
	}
	a <- (1 - level) / 2
	pct <- paste0(format(100 * c(a, 1 - a), trim = TRUE, digits = 3), " %")

	if (identical(object$type, "dyad")) {
		# jackknife: normal interval from the jackknife se (no replicate dist)
		z <- stats::qnorm(1 - a)
		ci <- cbind(object$point_est - z * object$se,
					object$point_est + z * object$se)
	} else {
		valid <- apply(object$coefs, 1, function(r) !any(is.na(r)))
		boot_valid <- object$coefs[valid, , drop = FALSE]
		if (nrow(boot_valid) < 2L) {
			cli::cli_abort("Bootstrap intervals require at least 2 valid replicates.")
		}
		ci <- t(apply(boot_valid, 2, quantile, probs = c(a, 1 - a)))
	}
	rownames(ci) <- object$param_names
	colnames(ci) <- pct

	if (!is.null(parm)) {
		ci <- ci[parm, , drop = FALSE]
	}
	ci
}
