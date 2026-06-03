# fit full-bilinear bipartite models

# compute the bipartite linear predictor
eta_tab_bipartite <- function(tab, W, W_recv, X, Z, p, p2, q) {
	d <- dim(X)
	n1 <- d[1]; n2 <- d[2]; Tt <- d[3]
	theta <- if (q > 0) tab[seq_len(q)] else numeric(0)
	alpha <- c(1, if (p > 1) tab[q + seq_len(p - 1)] else numeric(0))
	beta  <- tab[q + (p - 1) + seq_len(p2)]

	A <- matrix(matrix(W, n1 * n1, p) %*% alpha, n1, n1)
	B <- matrix(matrix(W_recv, n2 * n2, p2) %*% beta, n2, n2)

	eta <- array(0, dim = c(n1, n2, Tt))
	if (q > 0) {
		Zf <- if (length(dim(Z)) == 3) array(Z, dim = c(n1, n2, 1, Tt)) else Z
		for (k in seq_len(q)) {
			eta <- eta + theta[k] * array(Zf[, , k, ], dim = c(n1, n2, Tt))
		}
	}
	for (t in seq_len(Tt)) eta[, , t] <- eta[, , t] + A %*% X[, , t] %*% t(B)
	eta
}

# fit one flattened glm step
.bip_glm <- function(y, design, offset, family, family_obj) {
	df <- data.frame(.y = y, design, check.names = FALSE)
	rhs <- paste(sprintf("`%s`", colnames(design)), collapse = " + ")
	if (!is.null(offset)) {
		df$.off <- offset
		form <- stats::as.formula(paste0(".y ~ -1 + ", rhs, " + offset(.off)"))
	} else {
		form <- stats::as.formula(paste0(".y ~ -1 + ", rhs))
	}
	if (family == "normal") {
		fit <- stats::lm(form, data = df)
		dev <- sum(stats::residuals(fit)^2)
	} else {
		fit <- stats::glm(form, data = df, family = family_obj)
		dev <- stats::deviance(fit)
	}
	co <- stats::coef(fit)
	co[is.na(co)] <- 0
	list(coef = co, deviance = dev, converged = glm_converged(fit))
}

# fit the full-bilinear bipartite model
.fit_als_bipartite <- function(Y, W, W_recv, X, Z, family,
							   max_iter = 50, tol = 1e-6, trace = FALSE,
							   beta_init = NULL) {
	d <- dim(Y)
	n1 <- d[1]; n2 <- d[2]; Tt <- d[3]
	p  <- dim(W)[3]
	p2 <- dim(W_recv)[3]
	q  <- if (is.null(Z)) 0L else if (length(dim(Z)) == 3) 1L else dim(Z)[3]
	N  <- n1 * n2 * Tt

	family_obj <- switch(family,
		normal = NULL,
		poisson = stats::poisson(),
		binomial = stats::binomial(),
		cli::cli_abort("Unsupported family: {.val {family}}."))

	y <- as.vector(Y)

	# build exogenous design columns
	z_design <- NULL
	if (q > 0) {
		Zf <- if (length(dim(Z)) == 3) array(Z, dim = c(n1, n2, 1, Tt)) else Z
		z_design <- vapply(seq_len(q),
			function(k) as.vector(array(Zf[, , k, ], dim = c(n1, n2, Tt))),
			numeric(N))
		colnames(z_design) <- if (q > 0 && !is.null(dimnames(Z)[[3]])) {
			dimnames(Z)[[3]]
		} else {
			paste0("Z", seq_len(q))
		}
	}

	# stack bilinear influence values over time
	infl_vec <- function(L, R) {
		out <- array(0, dim = c(n1, n2, Tt))
		for (t in seq_len(Tt)) out[, , t] <- L %*% X[, , t] %*% t(R)
		as.vector(out)
	}

	# initialize the scaled influence vectors
	alpha <- c(1, rep(0, p - 1))
	beta  <- if (is.null(beta_init)) c(1, rep(0, p2 - 1)) else beta_init
	theta <- numeric(q)

	dev_old <- Inf
	dev_new <- NA_real_
	converged <- FALSE
	iter <- 0L
	for (iter in seq_len(max_iter)) {
		# hold beta and update theta and alpha
		B <- matrix(matrix(W_recv, n2 * n2, p2) %*% beta, n2, n2)
		off <- infl_vec(matrix(W[, , 1], n1, n1), B)
		a_design <- NULL
		if (p > 1) {
			a_design <- vapply(2:p,
				function(r) infl_vec(matrix(W[, , r], n1, n1), B), numeric(N))
			colnames(a_design) <- paste0("alphaW", 2:p)
		}
		des_a <- if (q > 0 && !is.null(a_design)) cbind(z_design, a_design)
			else if (q > 0) z_design else a_design
		if (!is.null(des_a)) {
			fa <- .bip_glm(y, des_a, off, family, family_obj)
			if (q > 0) theta <- fa$coef[seq_len(q)]
			if (p > 1) alpha <- c(1, fa$coef[q + seq_len(p - 1)])
		}

		# hold alpha and update theta and beta
		A <- matrix(matrix(W, n1 * n1, p) %*% alpha, n1, n1)
		b_design <- vapply(seq_len(p2),
			function(s) infl_vec(A, matrix(W_recv[, , s], n2, n2)), numeric(N))
		colnames(b_design) <- paste0("betaWr", seq_len(p2))
		des_b <- if (q > 0) cbind(z_design, b_design) else b_design
		fb <- .bip_glm(y, des_b, NULL, family, family_obj)
		if (q > 0) theta <- fb$coef[seq_len(q)]
		beta <- fb$coef[q + seq_len(p2)]

		dev_new <- fb$deviance
		if (trace) {
			cli::cli_alert_info("bipartite iter {.val {iter}}: deviance = {.val {sprintf('%.4f', dev_new)}}")
		}
		if (is.finite(dev_new) && is.finite(dev_old) &&
			abs(dev_old - dev_new) / (abs(dev_old) + 0.1) < tol) {
			converged <- TRUE
			break
		}
		dev_old <- dev_new
	}

	tab <- c(theta, if (p > 1) alpha[-1] else numeric(0), beta)
	list(theta = theta, alpha = alpha, beta = beta, tab = tab,
		 deviance = dev_new, iterations = iter, converged = converged,
		 p = p, p2 = p2, q = q)
}

# try several beta starts and keep the lowest deviance
.fit_als_bipartite_restart <- function(Y, W, W_recv, X, Z, family,
									   max_iter = 50, tol = 1e-6, trace = FALSE,
									   n_restarts = 5L) {
	p2 <- dim(W_recv)[3]
	best <- NULL
	starts <- c(list(NULL), if (n_restarts > 1) lapply(seq_len(n_restarts - 1),
		function(i) stats::rnorm(p2, sd = 0.5)) else list())
	for (b0 in starts) {
		fit <- tryCatch(
			.fit_als_bipartite(Y, W, W_recv, X, Z, family,
							   max_iter = max_iter, tol = tol, trace = trace,
							   beta_init = b0),
			error = function(e) NULL)
		if (is.null(fit) || !is.finite(fit$deviance)) next
		if (is.null(best) || fit$deviance < best$deviance) best <- fit
	}
	if (is.null(best)) {
		cli::cli_abort("Bipartite ALS failed on all restarts.")
	}
	if (!best$converged) {
		cli::cli_warn("Bipartite ALS did not converge within {.val {max_iter}} iterations (best of {n_restarts} restart{?s}).")
	}
	best
}

# build a sir_fit object for the full-bilinear bipartite model.
.sir_bipartite <- function(Y, W, W_recv, X, Z, family,
						   max_iter = 50, tol = 1e-6, trace = FALSE,
						   n_restarts = 5L) {
	d <- dim(Y)
	n1 <- d[1]; n2 <- d[2]; Tt <- d[3]

	# store one covariate with an explicit covariate axis
	if (!is.null(Z) && length(dim(Z)) == 3) Z <- array(Z, dim = c(n1, n2, 1, Tt))

	fit <- .fit_als_bipartite_restart(Y, W, W_recv, X, Z, family,
									  max_iter = max_iter, tol = tol, trace = trace,
									  n_restarts = n_restarts)
	p <- fit$p; p2 <- fit$p2; q <- fit$q

	eta <- eta_tab_bipartite(fit$tab, W, W_recv, X, Z, p, p2, q)
	fitted_values <- switch(family,
		poisson = exp(eta),
		binomial = 1 / (1 + exp(-eta)),
		normal = eta)

	A <- matrix(matrix(W, n1 * n1, p) %*% fit$alpha, n1, n1)
	B <- matrix(matrix(W_recv, n2 * n2, p2) %*% fit$beta, n2, n2)
	# carry sender/receiver labels onto the influence matrices when available
	send_nm <- dimnames(Y)[[1]]; recv_nm <- dimnames(Y)[[2]]
	if (!is.null(send_nm)) dimnames(A) <- list(send_nm, send_nm)
	if (!is.null(recv_nm)) dimnames(B) <- list(recv_nm, recv_nm)

	# parameter names follow the rest of the package
	w_names  <- if (!is.null(dimnames(W)[[3]])) dimnames(W)[[3]] else paste0("W", seq_len(p))
	wr_names <- if (!is.null(dimnames(W_recv)[[3]])) dimnames(W_recv)[[3]] else paste0("Wr", seq_len(p2))
	theta_names <- if (q > 0 && !is.null(dimnames(Z)[[3]])) dimnames(Z)[[3]]
		else if (q > 0) paste0("Z", seq_len(q)) else NULL
	alpha_names <- if (p > 1) paste0("(alphaW) ", w_names[-1]) else NULL
	beta_names  <- paste0("(betaWr) ", wr_names)

	summ <- data.frame(coef = fit$tab)
	summ$se <- NA; summ$rse <- NA; summ$t_se <- NA; summ$t_rse <- NA
	rownames(summ) <- c(theta_names, alpha_names, beta_names)

	ok <- !is.na(Y)
	nobs <- sum(ok)
	fitted_values[!ok] <- NA

	# build residual arrays
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

	sigma2 <- NULL
	sigma2_mle <- NULL
	if (family == "normal") {
		df_resid <- nobs - length(fit$tab)
		if (df_resid <= 0) {
			cli::cli_abort("Normal full-bilinear bipartite fit has non-positive residual degrees of freedom ({.val {df_resid}}).")
		}
		rss <- sum(response_resid[ok]^2)
		sigma2 <- rss / df_resid
		sigma2_mle <- rss / nobs
		pearson_resid <- response_resid / sqrt(sigma2)
	}

	# log-likelihood for AIC/BIC
	ll <- switch(family,
		poisson = sum(stats::dpois(Y[ok], lambda = pmax(fitted_values[ok], 1e-10), log = TRUE)),
		binomial = {
			pr <- pmin(pmax(fitted_values[ok], 1e-10), 1 - 1e-10)
			sum(stats::dbinom(Y[ok], 1, pr, log = TRUE))
		},
		normal = sum(stats::dnorm(Y[ok], mean = fitted_values[ok], sd = sqrt(sigma2_mle), log = TRUE)))

	result <- list(
		summ = summ,
		A = A,
		B = B,
		ll = ll,
		family = family,
		method = "ALS",
		tab = fit$tab,
		theta = fit$theta,
		alpha = fit$alpha,
		beta = fit$beta,
		p = p,
		p2 = p2,
		q = q,
		m = n1,
		n1 = n1,
		n2 = n2,
		bipartite = TRUE,
		full_bilinear = TRUE,
		n_periods = Tt,
		nobs = nobs,
		fitted.values = fitted_values,
		residuals = list(response = response_resid,
						 pearson = pearson_resid,
						 deviance = deviance_resid),
		vcov = NULL,
		vcov_robust = NULL,
		Y = Y,
		W = W,
		W_recv = W_recv,
		X = X,
		Z = Z,
		fix_receiver = FALSE,
		symmetric = FALSE,
		kron_mode = FALSE,
		dynamic_W = FALSE,
		iterations = fit$iterations,
		history = NULL,
		convergence = fit$converged,
		sigma2 = sigma2,
		sigma2_mle = sigma2_mle
	)
	class(result) <- c("sir_fit", "sir")
	result
}
