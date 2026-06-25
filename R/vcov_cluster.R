# cluster-robust variance for sir fits.
#
# the sandwich V = B M B clusters the per-observation scores on the actor margin
# (each cell's score is assigned to both endpoint actors) to account for the
# dyadic dependence that the classical and hc0 ses ignore. B is the inverse
# observed information (object$vcov), M is the actor-sum meat with an hc1
# g/(g-1) factor (g = #actors), and the result is psd-repaired.

# per-observation score matrix reproducing the c++ kernel's d_eta.
.sir_obs_scores <- function(object) {
	Y <- object$Y; W <- object$W; X <- object$X; Z <- object$Z
	family <- object$family; tab <- object$tab
	fix_receiver <- isTRUE(object$fix_receiver)

	# static W is m x m x p; dynamic W is m x m x p x T, so a_t/b_t vary by period
	W_dynamic <- !is.null(W) && length(dim(W)) == 4

	# symmetric (A = B) fits use a dedicated score with d eta / d gamma_r =
	# W_r X A' + A X W_r' over the upper-triangle off-diagonal cells
	if (isTRUE(object$symmetric) && identical(object$operator, "symmetric")) {
		return(.sir_obs_scores_symmetric(object))
	}

	n1 <- dim(X)[1]; n2 <- dim(X)[2]; T_len <- dim(X)[3]
	p  <- if (is.null(W)) 0L else dim(W)[3]
	zdim <- if (!is.null(Z)) length(dim(Z)) else 0L
	q  <- if (zdim == 4) dim(Z)[3] else if (zdim == 3) 1L else 0L
		drop_diagonal <- (n1 == n2) && !isTRUE(object$bipartite)
	sig2 <- if (family == "normal") {
		if (!is.null(object$sigma2)) object$sigma2 else 1
	} else 1

	# unpack parameters in the stored ordering.
	theta <- if (q > 0) tab[seq_len(q)] else numeric(0)
	if (fix_receiver) {
		# tab = [theta, alpha_1:p]; B = I, no beta.
		alpha <- if (p > 0) tab[q + seq_len(p)] else numeric(0)
		P <- q + p
	} else {
		# tab = [theta, alpha_2:p, beta_1:p]; alpha_1 = 1 fixed.
		alpha <- c(1, if (p > 1) tab[q + seq_len(p - 1)] else numeric(0))
		beta  <- if (p > 0) tab[q + (p - 1) + seq_len(p)] else numeric(0)
		P <- q + max(p - 1, 0) + p
	}

	# W_r for covariate r at period t: static W reuses one slice every t, dynamic
	# W indexes the period
	wget <- function(r, t) if (W_dynamic) W[, , r, t] else W[, , r]
	# A_t = sum_r alpha_r W_r(t); B_t = sum_r beta_r W_r(t) (B = I when fix_receiver)
	build_AB <- function(t) {
		A <- matrix(0, n1, n1)
		if (p > 0) for (r in seq_len(p)) A <- A + alpha[r] * wget(r, t)
		if (fix_receiver) {
			B <- diag(n2)
		} else {
			B <- matrix(0, n2, n2)
			if (p > 0) for (r in seq_len(p)) B <- B + beta[r] * wget(r, t)
		}
		list(A = A, B = B)
	}
	# static A/B are constant, so build once; dynamic rebuilds inside the loop
	if (!W_dynamic) { AB <- build_AB(1L); A <- AB$A; B <- AB$B }

	zget <- function(k, t) if (zdim == 4) Z[, , k, t] else Z[, , t]

	score_chunks <- vector("list", T_len)
	send_chunks  <- vector("list", T_len)
	recv_chunks  <- vector("list", T_len)

	for (t in seq_len(T_len)) {
		if (W_dynamic) { AB <- build_AB(t); A <- AB$A; B <- AB$B }
		Xt <- X[, , t]; Yt <- Y[, , t]
		eta <- matrix(0, n1, n2)
		if (q > 0) for (k in seq_len(q)) eta <- eta + theta[k] * zget(k, t)
		AX <- A %*% Xt
		eta <- eta + AX %*% t(B)
		mu <- switch(family,
			poisson  = exp(eta),
			binomial = pmin(pmax(1 / (1 + exp(-eta)), 1e-10), 1 - 1e-10),
			normal   = eta)
		resid <- (Yt - mu) / sig2

		# d_eta columns, in the stored parameter order; each is an n1 x n2 matrix.
		deta <- vector("list", P); cc <- 1L
		if (q > 0) for (k in seq_len(q)) { deta[[cc]] <- zget(k, t); cc <- cc + 1L }
		if (fix_receiver) {
			# alpha_r (r = 1..p): d eta / d alpha_r = W_r X_t   (B = I)
			if (p > 0) for (r in seq_len(p)) { deta[[cc]] <- wget(r, t) %*% Xt; cc <- cc + 1L }
		} else {
			# alpha_r (r = 2..p): d eta / d alpha_r = W_r X_t B'
			if (p > 1) for (r in 2:p) { deta[[cc]] <- wget(r, t) %*% Xt %*% t(B); cc <- cc + 1L }
			# beta_r (r = 1..p): d eta / d beta_r = A X_t W_r'
			if (p > 0) for (r in seq_len(p)) { deta[[cc]] <- AX %*% t(wget(r, t)); cc <- cc + 1L }
		}

		rvec <- as.vector(resid)
		sc <- matrix(0, n1 * n2, P)
		for (cidx in seq_len(P)) sc[, cidx] <- rvec * as.vector(deta[[cidx]])

		ij <- expand.grid(i = seq_len(n1), j = seq_len(n2))
		score_chunks[[t]] <- sc
		send_chunks[[t]]  <- ij$i
		recv_chunks[[t]]  <- ij$j
	}

	scores <- do.call(rbind, score_chunks)
	sender <- unlist(send_chunks)
	receiver <- unlist(recv_chunks)

		# keep only modeled observations: off-diagonal for one-mode square fits and
		# finite y. Square bipartite fits keep diagonal sender-receiver cells.
		yvec <- unlist(lapply(seq_len(T_len), function(t) as.vector(Y[, , t])))
		keep <- is.finite(yvec) & rowSums(!is.finite(scores)) == 0
		if (drop_diagonal) keep <- keep & (sender != receiver)

	list(scores = scores[keep, , drop = FALSE],
		 sender = sender[keep], receiver = receiver[keep])
}

# per-observation scores for the symmetric (A = B) operator over the upper-
# triangle off-diagonal cells. returns scores plus the actor
# index for each endpoint (i, j) and the time index.
.sir_obs_scores_symmetric <- function(object) {
	Y <- object$Y; W <- object$W; X <- object$X; Z <- object$Z
	family <- object$family; tab <- object$tab
	n <- dim(X)[1]; T_len <- dim(X)[3]
	p <- if (is.null(W)) 0L else dim(W)[3]
	zdim <- if (!is.null(Z)) length(dim(Z)) else 0L
	q  <- if (zdim == 4) dim(Z)[3] else if (zdim == 3) 1L else 0L
	sig2 <- if (family == "normal" && !is.null(object$sigma2)) object$sigma2 else 1
	theta <- if (q > 0) tab[seq_len(q)] else numeric(0)
	gamma <- tab[q + seq_len(p)]
	W_dynamic <- length(dim(W)) == 4
	wget <- function(r, t) if (W_dynamic) W[, , r, t] else W[, , r]
	# A_t = sum_r gamma_r W_r(t); constant across t for static W
	build_A <- function(t) {
		Wt <- if (W_dynamic) W[, , , t] else W
		matrix(matrix(Wt, n * n, p) %*% gamma, n, n)
	}
	if (!W_dynamic) A <- build_A(1L)
	zget <- function(k, t) if (zdim == 4) Z[, , k, t] else Z[, , t]
	P <- q + p
	um <- upper.tri(matrix(0, n, n))
	ut <- which(um)
	ij <- expand.grid(i = seq_len(n), j = seq_len(n))[ut, ]

	chunks <- vector("list", T_len)
	keepL <- vector("list", T_len)
	for (t in seq_len(T_len)) {
		if (W_dynamic) A <- build_A(t)
		Xt <- X[, , t]; Yt <- Y[, , t]
		eta <- A %*% Xt %*% t(A)
		if (q > 0) for (k in seq_len(q)) eta <- eta + theta[k] * zget(k, t)
		mu <- switch(family,
			poisson = exp(pmin(pmax(eta, -500), 500)),
			binomial = pmin(pmax(1 / (1 + exp(-eta)), 1e-10), 1 - 1e-10),
			normal = eta)
		resid <- (Yt - mu) / sig2
		deta <- vector("list", P); cc <- 1L
		if (q > 0) for (k in seq_len(q)) { deta[[cc]] <- zget(k, t); cc <- cc + 1L }
		AX <- A %*% Xt; XtA <- Xt %*% t(A)
		for (r in seq_len(p)) {
			Wr <- wget(r, t)
			deta[[cc]] <- Wr %*% XtA + AX %*% t(Wr); cc <- cc + 1L
		}
		sc <- matrix(0, length(ut), P)
		rvec <- resid[ut]
		for (cidx in seq_len(P)) sc[, cidx] <- rvec * deta[[cidx]][ut]
		ok <- is.finite(Yt[ut]) & rowSums(!is.finite(sc)) == 0
		chunks[[t]] <- sc[ok, , drop = FALSE]
		keepL[[t]] <- cbind(i = ij$i[ok], j = ij$j[ok])
	}
	scores <- do.call(rbind, chunks)
	idx <- do.call(rbind, keepL)
	list(scores = scores, sender = idx[, "i"], receiver = idx[, "j"], symmetric = TRUE)
}

# meat for one clustering: crossprod of the per-group score sums.
.cluster_meat <- function(scores, group) {
	totals <- rowsum(scores, as.integer(factor(group)))
	crossprod(totals)
}

# PSD repair: symmetrize and floor eigenvalues at zero.
.psd_repair <- function(V) {
	V <- (V + t(V)) / 2
	ev <- eigen(V, symmetric = TRUE)
	lam <- pmax(ev$values, 0)
	ev$vectors %*% (lam * t(ev$vectors))
}

# the cluster-robust covariance: assign each cell's score to both endpoint
# actors, sandwich with the classical bread, and return a P x P covariance in
# the same ordering as object$vcov.
.sir_vcov_cluster <- function(object) {
	bread <- object$vcov
	if (is.null(bread)) {
		cli::cli_abort(c(
			"Cluster-robust SEs require the classical covariance (the bread).",
			"i" = "Refit with {.code calc_se = TRUE}."
		))
	}
	if (isFALSE(object$se_reliable)) {
		cli::cli_warn(c(
			"Cluster-robust SEs use an unstable Hessian as the bread for this fit.",
			"i" = "Treat them as sensitivity checks; consider refitting, simplifying the model, or using {.code boot_sir(type = \"dyad\")}."
		))
	}
	sc <- .sir_obs_scores(object)
	scores <- sc$scores
	if (nrow(scores) == 0L) cli::cli_abort("No modeled observations available for clustering.")
	# guard a dimension mismatch (e.g. an unexpected parameter layout).
	if (ncol(scores) != nrow(bread)) {
		cli::cli_abort(c(
			"Score dimension ({ncol(scores)}) does not match the covariance ({nrow(bread)}).",
			"i" = "Cluster-robust SEs are unavailable for this fit; use {.code boot_sir()}."
		))
	}

	# actor-margin cluster-robust sandwich, used for both directed and symmetric
	# fits: assign every cell's score to BOTH endpoint actors and sum within
	# actor, meat = sum_a g_a g_a' with an HC1 G/(G-1) factor (G = #actors). the
	# influence-coefficient score loads on the actor margin, so this is the
	# covariance that matters.
	sc2 <- rbind(scores, scores)
	# in a bipartite fit senders and receivers are distinct populations, so tag
	# the two margins apart (sender i and receiver i are different actors); in a
	# one-mode/symmetric fit node i is the same actor on both sides
	actor <- if (isTRUE(object$bipartite))
		c(paste0("s", sc$sender), paste0("r", sc$receiver))
	else c(sc$sender, sc$receiver)
	meat <- .cluster_meat(sc2, actor)
	G <- length(unique(actor))
	if (G > 1) meat <- meat * (G / (G - 1))
	V <- bread %*% meat %*% bread
	V <- .psd_repair(V)
	dn <- dimnames(bread)
	if (!is.null(dn)) dimnames(V) <- dn
	# carry the cluster count so confint can use a t(G-1) reference
	attr(V, "cluster_df") <- G - 1L
	V
}
