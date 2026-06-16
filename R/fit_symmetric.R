# undirected fit: A = B = sum_k gamma_k W_k, so the bilinear term is the
# quadratic form A X A', symmetric by construction. all p gammas are free (the
# quadratic form's scale is data-identified); A is identified only up to global
# sign, fixed so the largest-magnitude gamma is positive. fits by BFGS on the
# upper-triangle off-diagonal cells.

# operator A_t = sum_k gamma_k W_{k,t}. static W (3D) gives the same A every
# period; dynamic W (4D, n x n x p x T) gives a per-period A.
.sym_A_t <- function(W, gamma, n, p, t, dynamic) {
	if (dynamic) matrix(matrix(W[, , , t], n * n, p) %*% gamma, n, n)
	else matrix(matrix(W, n * n, p) %*% gamma, n, n)
}

# linear predictor. tab packs [theta, gamma_1..gamma_p] (all gammas free).
# W slices must have a zero diagonal. W may be static (3D) or dynamic (4D).
eta_tab_symmetric <- function(tab, W, X, Z, p, q) {
	d <- dim(X)
	n <- d[1]; Tt <- d[3]
	dynamic <- length(dim(W)) == 4
	theta <- if (q > 0) tab[seq_len(q)] else numeric(0)
	gamma <- tab[q + seq_len(p)]
	eta <- array(0, dim = c(n, n, Tt))
	if (q > 0) {
		Zf <- if (length(dim(Z)) == 3) array(Z, dim = c(n, n, 1, Tt)) else Z
		for (k in seq_len(q)) eta <- eta + theta[k] * array(Zf[, , k, ], dim = c(n, n, Tt))
	}
	A <- if (!dynamic) .sym_A_t(W, gamma, n, p, 1L, FALSE) else NULL
	for (t in seq_len(Tt)) {
		At <- if (dynamic) .sym_A_t(W, gamma, n, p, t, TRUE) else A
		eta[, , t] <- eta[, , t] + At %*% X[, , t] %*% t(At)
	}
	eta
}

# one BFGS fit from a given start. returns gamma, theta, deviance, convergence,
# gradient norm. par = [theta, gamma_1..p].
.fit_symmetric_once <- function(Y, W, X, Z, family, p, q, n, Tt,
								par0, max_iter = 200, tol = 1e-8) {
	dynamic <- length(dim(W)) == 4
	um <- upper.tri(matrix(0, n, n))                 # one cell per unordered dyad
	# slice covariate k at time t (works for static or dynamic W)
	Wkt <- function(k, t) if (dynamic) W[, , k, t] else W[, , k]
	# z design slices (n x n x Tt) per covariate
	Zf <- if (q > 0) {
		if (length(dim(Z)) == 3) array(Z, dim = c(n, n, 1, Tt)) else Z
	} else NULL

	gamma_of <- function(par) par[q + seq_len(p)]
	eta_of <- function(par) {
		theta <- if (q > 0) par[seq_len(q)] else numeric(0)
		gamma <- gamma_of(par)
		eta <- array(0, dim = c(n, n, Tt))
		if (q > 0) for (k in seq_len(q)) eta <- eta + theta[k] * array(Zf[, , k, ], dim = c(n, n, Tt))
		A <- if (!dynamic) .sym_A_t(W, gamma, n, p, 1L, FALSE) else NULL
		for (t in seq_len(Tt)) {
			At <- if (dynamic) .sym_A_t(W, gamma, n, p, t, TRUE) else A
			eta[, , t] <- eta[, , t] + At %*% X[, , t] %*% t(At)
		}
		eta
	}

	# negative log-likelihood over upper-triangle off-diagonal non-missing cells
	nll <- function(par) {
		eta <- eta_of(par)
		v <- 0
		for (t in seq_len(Tt)) {
			yt <- Y[, , t][um]; et <- eta[, , t][um]
			ok <- is.finite(yt) & is.finite(et)
			yt <- yt[ok]; et <- et[ok]
			if (!length(yt)) next
			if (family == "normal") {
				v <- v + 0.5 * sum((yt - et)^2)
			} else if (family == "poisson") {
				et <- pmin(pmax(et, -500), 500)        # clamp before exp
				v <- v + sum(exp(et) - yt * et)
			} else {                                    # binomial
				et <- pmin(pmax(et, -500), 500)
				v <- v + sum(log1p(exp(et)) - yt * et)  # -loglik of bernoulli logit
			}
		}
		if (!is.finite(v)) return(1e10)                # finite penalty for restarts
		v
	}

	# analytic gradient. d eta / d gamma_r = W_r X A' + A X W_r' (r >= 2);
	# d eta / d theta_k = Z_k; working residual is mu - y for the family.
	grad <- function(par) {
		theta <- if (q > 0) par[seq_len(q)] else numeric(0)
		gamma <- gamma_of(par)
		g <- numeric(length(par))
		A <- if (!dynamic) .sym_A_t(W, gamma, n, p, 1L, FALSE) else NULL
		# the derivative operators depend on X (and on A_t for dynamic W)
		for (t in seq_len(Tt)) {
			Xt <- X[, , t]
			At <- if (dynamic) .sym_A_t(W, gamma, n, p, t, TRUE) else A
			eta_t <- At %*% Xt %*% t(At)
			if (q > 0) for (k in seq_len(q)) eta_t <- eta_t + theta[k] * Zf[, , k, t]
			yt <- Y[, , t]
			w_t <- switch(family,
				normal   = eta_t - yt,
				poisson  = exp(pmin(pmax(eta_t, -500), 500)) - yt,
				binomial = 1 / (1 + exp(-pmin(pmax(eta_t, -500), 500))) - yt)
			# keep only upper-triangle off-diagonal cells; zero the rest
			mask <- um & is.finite(w_t) & is.finite(eta_t)
			w_use <- matrix(0, n, n); w_use[mask] <- w_t[mask]
			if (q > 0) for (k in seq_len(q)) {
				g[k] <- g[k] + sum(w_use * Zf[, , k, t])
			}
			AX <- At %*% Xt
			XtA <- Xt %*% t(At)
			for (r in seq_len(p)) {
				Wr <- Wkt(r, t)
				d_eta_r <- Wr %*% XtA + AX %*% t(Wr)   # W_r X A' + A X W_r'
				g[q + r] <- g[q + r] + sum(w_use * d_eta_r)
			}
		}
		g
	}

	fit <- tryCatch(
		stats::optim(par0, nll, grad, method = "BFGS",
					 control = list(maxit = max_iter, reltol = tol)),
		error = function(e) NULL)
	if (is.null(fit)) return(NULL)
	gn <- sqrt(sum(grad(fit$par)^2))
	list(par = fit$par, deviance = 2 * fit$value, convergence = fit$convergence,
		 grad_norm = gn,
		 gamma = fit$par[q + seq_len(p)],
		 theta = if (q > 0) fit$par[seq_len(q)] else numeric(0))
}

# one-step gauss-newton start for gamma_2..p (gamma_1 starts at 1): regress the
# residual at A = W[,,1]
# on the gradient-direction basis.
.symmetric_init <- function(Y, W, X, Z, family, p, q, n, Tt) {
	# warm start on the linear-predictor scale
	dynamic <- length(dim(W)) == 4
	Wkt <- function(k, t) if (dynamic) W[, , k, t] else W[, , k]
	um <- which(upper.tri(matrix(0, n, n)))
	# link-scale target minus the anchor prediction A1 X A1' (A1 from W[,,1])
	rows <- list(); rhs <- c()
	for (t in seq_len(Tt)) {
		Xt <- X[, , t]
		A1 <- Wkt(1, t)
		base <- A1 %*% Xt %*% t(A1)
		yt <- Y[, , t]
		link_y <- switch(family,
			normal = yt,
			poisson = log(pmax(yt, 0) + 1),
			binomial = stats::qlogis(pmin(pmax(yt, 0.01), 0.99)))
		resid <- (link_y - base)[um]
		Dt <- matrix(0, length(um), q + max(p - 1, 0))
		if (q > 0) {
			Zf <- if (length(dim(Z)) == 3) array(Z, dim = c(n, n, 1, Tt)) else Z
			for (k in seq_len(q)) Dt[, k] <- Zf[, , k, t][um]
		}
		if (p > 1) {
			AX <- A1 %*% Xt; XtA <- Xt %*% t(A1)
			for (r in 2:p) {
				Wr <- Wkt(r, t)
				Dt[, q + r - 1] <- (Wr %*% XtA + AX %*% t(Wr))[um]
			}
		}
		ok <- is.finite(resid) & apply(Dt, 1, function(z) all(is.finite(z)))
		rhs <- c(rhs, resid[ok]); rows[[t]] <- Dt[ok, , drop = FALSE]
	}
	# assemble the full start [theta, gamma_1 = 1, gamma_2..p]; the gauss-newton
	# step estimates gamma_2..p relative to the anchor, gamma_1 starts at 1
	assemble <- function(v) c(if (q > 0) v[seq_len(q)] else numeric(0), 1,
							  if (p > 1) v[q + seq_len(p - 1)] else numeric(0))
	D <- do.call(rbind, rows)
	if (is.null(D) || nrow(D) < ncol(D)) return(assemble(numeric(q + max(p - 1, 0))))
	coef <- tryCatch(as.numeric(stats::lsfit(D, rhs, intercept = FALSE)$coefficients),
					 error = function(e) numeric(ncol(D)))
	coef[!is.finite(coef)] <- 0
	assemble(coef)
}

# run a warm start plus jittered restarts and keep the lowest-deviance fit
.fit_symmetric <- function(Y, W, X, Z, family, max_iter = 200, tol = 1e-8,
						   n_restarts = 5L, trace = FALSE) {
	d <- dim(Y); n <- d[1]; Tt <- d[3]
	p <- dim(W)[3]
	q <- if (is.null(Z)) 0L else if (length(dim(Z)) == 3) 1L else dim(Z)[3]
	npar <- q + p

	starts <- list()
	starts[[1]] <- .symmetric_init(Y, W, X, Z, family, p, q, n, Tt)
	if (n_restarts > 1 && npar > 0) {
		for (i in 2:n_restarts) {
			# small jitter; large jumps diverge in the quadratic form
			starts[[i]] <- starts[[1]] + stats::rnorm(npar, sd = 0.05)
		}
	} else if (npar == 0) {
		starts <- list(numeric(0))
	}

	best <- NULL
	spread <- c()
	for (s in starts) {
		f <- .fit_symmetric_once(Y, W, X, Z, family, p, q, n, Tt,
								 par0 = s, max_iter = max_iter, tol = tol)
		if (is.null(f) || !is.finite(f$deviance)) next
		spread <- c(spread, f$deviance)
		if (is.null(best) || f$deviance < best$deviance) best <- f
		if (trace) cli::cli_alert_info("symmetric restart deviance = {.val {sprintf('%.4f', f$deviance)}}")
	}
	if (is.null(best)) cli::cli_abort("Symmetric fit failed on all restarts.")

	# A X A' is invariant to A -> -A, so fix the global sign so the largest-
	# magnitude gamma is positive (stable even when gamma_1 is near zero); the
	# flip leaves the deviance and fitted values unchanged
	g_idx <- q + seq_len(p)
	anchor <- which.max(abs(best$par[g_idx]))
	if (best$par[g_idx][anchor] < 0) {
		best$par[g_idx] <- -best$par[g_idx]
		best$gamma <- -best$gamma
	}

	# converged if optim succeeded and the gradient is small relative to the
	# objective scale, or all restarts agree on the deviance
	grad_tol <- 1e-4 * max(abs(best$deviance), 1)
	restarts_agree <- length(spread) >= 2 &&
		all(abs(spread[is.finite(spread)] - best$deviance) /
			(abs(best$deviance) + 0.1) < 1e-3)
	best$converged <- isTRUE(best$convergence == 0) &&
		(best$grad_norm < grad_tol || restarts_agree)
	# flag competing optima: a restart converged to a materially worse deviance
	# than the best one (best is the minimum, so look for higher deviances)
	fin <- spread[is.finite(spread)]
	multimodal <- length(fin) >= 2 &&
		(max(fin) - best$deviance) / (abs(best$deviance) + 0.1) > 1e-3
	if (!best$converged) {
		cli::cli_warn(c(
			"Symmetric fit did not converge cleanly (grad norm {.val {sprintf('%.2e', best$grad_norm)}}).",
			"i" = "Try more {.arg n_restarts} or check that {.arg W} is well-conditioned."
		))
	} else if (multimodal) {
		cli::cli_warn("Symmetric fit is multimodal: restarts found competing optima. Inspect {.code fit$restart_spread}.")
	}
	best$tab <- best$par
	best$p <- p; best$q <- q; best$n <- n; best$Tt <- Tt
	best$restart_spread <- spread
	best
}

# build a sir_fit object for the symmetric (A = B) model.
.sir_symmetric <- function(Y, W, X, Z, family, max_iter = 200, tol = 1e-8,
						   n_restarts = 5L, trace = FALSE, calc_se = TRUE) {
	d <- dim(Y); n <- d[1]; Tt <- d[3]

	# fit each unordered dyad once: blank the lower triangle and the diagonal
	for (t in seq_len(Tt)) {
		Yt <- Y[, , t]
		Yt[lower.tri(Yt)] <- NA
		diag(Yt) <- NA
		Y[, , t] <- Yt
	}
	if (!is.null(Z) && length(dim(Z)) == 3) Z <- array(Z, dim = c(n, n, 1, Tt))

	fit <- .fit_symmetric(Y, W, X, Z, family, max_iter = max_iter, tol = tol,
						  n_restarts = n_restarts, trace = trace)
	p <- fit$p; q <- fit$q
	dynamic <- length(dim(W)) == 4

	# static W gives one operator A; dynamic W gives a per-period A_t array
	nm <- dimnames(Y)[[1]]
	if (dynamic) {
		A <- array(0, dim = c(n, n, Tt))
		for (t in seq_len(Tt)) A[, , t] <- .sym_A_t(W, fit$gamma, n, p, t, TRUE)
		if (!is.null(nm)) dimnames(A) <- list(nm, nm, NULL)
	} else {
		A <- .sym_A_t(W, fit$gamma, n, p, 1L, FALSE)
		if (!is.null(nm)) dimnames(A) <- list(nm, nm)
	}

	# fitted values: full symmetric matrix on the response scale, NA diagonal
	eta <- eta_tab_symmetric(fit$tab, W, X, Z, p, q)
	fitted_values <- switch(family,
		poisson = exp(pmin(pmax(eta, -500), 500)),
		binomial = 1 / (1 + exp(-pmin(pmax(eta, -500), 500))),
		normal = eta)
	for (t in seq_len(Tt)) diag(fitted_values[, , t]) <- NA

	# all p gammas are estimated and reported; direct effects carry the (Z) tag
	# for parity with the directed fit
	w_names <- if (!is.null(dimnames(W)[[3]])) dimnames(W)[[3]] else paste0("W", seq_len(p))
	theta_names <- if (q > 0 && !is.null(dimnames(Z)[[3]])) paste0("(Z) ", dimnames(Z)[[3]])
		else if (q > 0) paste0("(Z) Z", seq_len(q)) else NULL
	gamma_names <- paste0("(gammaW) ", w_names)

	summ <- data.frame(coef = fit$tab)
	summ$se <- NA; summ$rse <- NA; summ$t_se <- NA; summ$t_rse <- NA
	rownames(summ) <- c(theta_names, gamma_names)

	# count on the upper-triangle off-diagonal cells actually used
	ok <- !is.na(Y) & !is.na(fitted_values)
	nobs <- sum(ok)

	response_resid <- Y - fitted_values
	pearson_resid <- switch(family,
		poisson = response_resid / sqrt(fitted_values),
		binomial = response_resid / sqrt(fitted_values * (1 - fitted_values)),
		normal = response_resid)
	deviance_resid <- response_resid
	if (family == "poisson") {
		dc <- 2 * (ifelse(Y > 0, Y * log(Y / fitted_values), 0) - (Y - fitted_values))
		dc[dc < 0] <- 0
		deviance_resid <- sign(response_resid) * sqrt(dc)
	} else if (family == "binomial") {
		dc <- 2 * (ifelse(Y > 0, Y * log(Y / fitted_values), 0) +
				   ifelse(Y < 1, (1 - Y) * log((1 - Y) / (1 - fitted_values)), 0))
		dc[dc < 0] <- 0
		deviance_resid <- sign(response_resid) * sqrt(dc)
	}

	sigma2 <- NULL; sigma2_mle <- NULL
	if (family == "normal") {
		df_resid <- nobs - length(fit$tab)
		if (df_resid <= 0) cli::cli_abort("Symmetric normal fit has non-positive residual df ({.val {df_resid}}).")
		rss <- sum(response_resid[ok]^2)
		sigma2 <- rss / df_resid
		sigma2_mle <- rss / nobs
		pearson_resid <- response_resid / sqrt(sigma2)
	}

	ll <- switch(family,
		poisson = sum(stats::dpois(Y[ok], lambda = pmax(fitted_values[ok], 1e-10), log = TRUE)),
		binomial = {
			pr <- pmin(pmax(fitted_values[ok], 1e-10), 1 - 1e-10)
			sum(stats::dbinom(Y[ok], 1, pr, log = TRUE))
		},
		normal = sum(stats::dnorm(Y[ok], mean = fitted_values[ok], sd = sqrt(sigma2_mle), log = TRUE)))

	# spectral gain rho(A)^2/(n-1) governs the recursion's stationarity; for
	# dynamic W take the largest spectral radius across periods. only the
	# poisson/normal recursion can explode -- the binomial outcome is bounded, so
	# a large gain there is not a stationarity problem
	rho_one <- function(M) tryCatch(max(abs(eigen(M, only.values = TRUE)$values)), error = function(e) NA_real_)
	rho_A <- if (dynamic) max(vapply(seq_len(Tt), function(t) rho_one(A[, , t]), numeric(1))) else rho_one(A)
	gain <- if (is.finite(rho_A)) rho_A^2 / max(n - 1, 1) else NA_real_
	gain_matters <- family %in% c("poisson", "normal")
	if (gain_matters && is.finite(gain) && gain >= 1) {
		cli::cli_warn("Estimated operator is non-stationary: spectral gain rho(A)^2/(n-1) = {.val {sprintf('%.2f', gain)}} >= 1.")
	}

	# classical covariance: inverse Hessian of the NLL at the optimum, also the
	# bread for the cluster sandwich
	vcov_mat <- NULL; se_reliable <- NA
	if (isTRUE(calc_se) && length(fit$tab) > 0) {
		vcov_mat <- tryCatch({
			H <- numDeriv::hessian(.sym_nll_for_hess, fit$tab,
								   Y = Y, W = W, X = X, Z = Z, family = family,
								   p = p, q = q, n = n, Tt = Tt, sigma2 = sigma2)
			Hs <- (H + t(H)) / 2
			ev <- eigen(Hs, symmetric = TRUE, only.values = TRUE)$values
			se_reliable <- all(ev > 0) && (max(ev) / min(ev) < 1e10)
			V <- tryCatch(solve(Hs), error = function(e) MASS::ginv(Hs))
			V <- (V + t(V)) / 2
			dimnames(V) <- list(rownames(summ), rownames(summ))
			V
		}, error = function(e) NULL)
		if (!is.null(vcov_mat)) {
			dv <- diag(vcov_mat)
			se <- sqrt(ifelse(dv >= 0, dv, NA_real_))
			summ$se <- se
			summ$t_se <- summ$coef / se
		}
	}

	result <- list(
		summ = summ, A = A, B = A, ll = ll, family = family, method = "optim",
		tab = fit$tab, theta = fit$theta, alpha = fit$gamma, beta = fit$gamma,
		gamma = fit$gamma, p = p, q = q, m = n, n1 = n, n2 = n,
		bipartite = FALSE, symmetric = TRUE, full_bilinear = FALSE,
		operator = "symmetric", rho_A = rho_A, gain = gain,
		stationary = if (gain_matters) !(is.finite(gain) && gain >= 1) else NA,
		n_periods = Tt, nobs = nobs, fitted.values = fitted_values,
		residuals = list(response = response_resid, pearson = pearson_resid,
						 deviance = deviance_resid),
		vcov = vcov_mat, vcov_robust = NULL, se_reliable = se_reliable,
		Y = Y, W = W, X = X, Z = Z,
		fix_receiver = FALSE, kron_mode = FALSE, dynamic_W = dynamic,
		iterations = NA, history = NULL, convergence = fit$converged,
		restart_spread = fit$restart_spread, sigma2 = sigma2)
	class(result) <- c("sir_fit", "sir")
	result
}

# upper-triangle off-diagonal NLL in tab = [theta, gamma_1..p], used for the
# observed-information Hessian. normal scales by sigma2; others use the GLM NLL.
.sym_nll_for_hess <- function(tab, Y, W, X, Z, family, p, q, n, Tt, sigma2) {
	eta <- eta_tab_symmetric(tab, W, X, Z, p, q)
	um <- upper.tri(matrix(0, n, n))
	v <- 0
	for (t in seq_len(Tt)) {
		yt <- Y[, , t][um]; et <- eta[, , t][um]
		ok <- is.finite(yt) & is.finite(et)
		yt <- yt[ok]; et <- et[ok]
		if (!length(yt)) next
		if (family == "normal") {
			s2 <- if (is.null(sigma2) || !is.finite(sigma2) || sigma2 <= 0) 1 else sigma2
			v <- v + sum(0.5 * (yt - et)^2 / s2)
		} else if (family == "poisson") {
			et <- pmin(pmax(et, -500), 500)
			v <- v + sum(exp(et) - yt * et)
		} else {
			et <- pmin(pmax(et, -500), 500)
			v <- v + sum(log1p(exp(et)) - yt * et)
		}
	}
	v
}
