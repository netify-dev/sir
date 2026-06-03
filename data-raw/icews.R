# derive bundled icews example data

raw_path = Sys.getenv("SIR_ICEWS_RAW")
if (!nzchar(raw_path)) {
	stop(
		"set SIR_ICEWS_RAW to the path for replArchive/data/socRegData.rda",
		call. = FALSE
	)
}
if (!file.exists(raw_path)) {
	stop("SIR_ICEWS_RAW does not point to an existing file", call. = FALSE)
}

data_env = new.env()
load(raw_path, envir = data_env)
required_objects = c("Y", "W_avg", "X", "Z")
missing_objects = setdiff(required_objects, ls(data_env))
if (length(missing_objects)) {
	stop(
		"socRegData.rda is missing required object(s): ",
		paste(missing_objects, collapse = ", "),
		call. = FALSE
	)
}
Y = data_env$Y
W_avg = data_env$W_avg
X = data_env$X
Z = data_env$Z
if (length(dim(Y)) != 3 || length(dim(X)) != 3 || length(dim(W_avg)) != 3 ||
	length(dim(Z)) != 4) {
	stop("unexpected source dimensions for Y, X, W_avg, or Z", call. = FALSE)
}
if (!all(dim(Y) == dim(X))) {
	stop("Y and X must have matching dimensions", call. = FALSE)
}

# rank countries by conflict involvement
total_conflict = apply(Y, c(1, 2), sum, na.rm = TRUE)
activity = rowSums(total_conflict) + colSums(total_conflict)
keep = sort(order(activity, decreasing = TRUE)[seq_len(50)])

# subset outcome and lagged state
Y_subset = Y[keep, keep, , drop = FALSE]
X_subset = X[keep, keep, , drop = FALSE]

# build influence covariates
w_keep = c("int", "ally", "verbCoop", "minDistLog")
if (!all(w_keep %in% dimnames(W_avg)[[3]])) {
	stop("W_avg is missing required covariate slice(s)", call. = FALSE)
}
W_subset = W_avg[keep, keep, w_keep, drop = FALSE]
W_subset[, , "verbCoop"] = log(W_subset[, , "verbCoop"] + 1)

# build direct-effect covariates
z_keep = c("mConf", "mConf_ji", "minDistLog", "ally", "verbCoop")
if (!all(z_keep %in% dimnames(Z)[[3]])) {
	stop("Z is missing required covariate slice(s)", call. = FALSE)
}
Z_subset = Z[keep, keep, z_keep, , drop = FALSE]
Z_subset[, , "verbCoop", ] = log(Z_subset[, , "verbCoop", ] + 1)

icews = list(
	Y = Y_subset,
	X = X_subset,
	W = W_subset,
	Z = Z_subset,
	countries = dimnames(Y_subset)[[1]],
	dates = dimnames(Y_subset)[[3]],
	metadata = list(
		source = list(
			paper_doi = "10.1017/pan.2025.10013",
			replication_doi = "10.7910/DVN/VTFDX6",
			raw_file = "replArchive/data/socRegData.rda",
			repository = "https://github.com/s7minhas/sir_paper",
			access_date = "2026-06-02"
		),
		bundled = list(
			country_rule = "Top 50 countries by total sent plus received material-conflict counts, sorted into source order",
			dates = c(first = dimnames(Y_subset)[[3]][1], last = tail(dimnames(Y_subset)[[3]], 1)),
			transforms = c(
				"X is log(Y_{t-1} + 1), with the first bundled period using the retained January 2005 lag",
				"W and Z verbCoop are log(x + 1)",
				"one-mode diagonals are stored as zero and ignored by sir()"
			)
		)
	)
)

stopifnot(
	identical(dim(icews$Y), c(50L, 50L, 95L)),
	identical(dim(icews$X), c(50L, 50L, 95L)),
	identical(dim(icews$W), c(50L, 50L, 4L)),
	identical(dim(icews$Z), c(50L, 50L, 5L, 95L)),
	identical(icews$dates[1], "2005-02-01"),
	identical(tail(icews$dates, 1), "2012-12-01")
)

if (requireNamespace("usethis", quietly = TRUE)) {
	usethis::use_data(icews, overwrite = TRUE, compress = "xz")
} else {
	save(icews, file = file.path("data", "icews.rda"), compress = "xz")
}
