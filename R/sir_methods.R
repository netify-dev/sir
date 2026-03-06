#' @importFrom cli cli_h1 cli_h2 cli_text cli_ul cli_alert_success cli_alert_warning cli_rule cli_abort cli_warn cli_inform
#' @importFrom stats printCoefmat deviance quantile pnorm symnum AIC BIC logLik fitted residuals coef confint nobs
#' @keywords internal
NULL

#' Extract Model Coefficients from a SIR Model
#'
#' Returns the estimated parameter values from a fitted SIR model. These
#' include exogenous covariate effects (theta), sender influence weights
#' (alpha), and receiver influence weights (beta).
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Additional arguments (unused).
#' @return Named numeric vector of estimated coefficients.
#' @seealso \code{\link{confint.sir}} for confidence intervals,
#'   \code{\link{vcov.sir}} for the variance-covariance matrix.
#' @export
coef.sir <- function(object, ...) {
  object$summ$coef
}

#' Extract Fitted Values from a SIR Model
#'
#' Returns the fitted values on the response scale: expected counts for Poisson,
#' probabilities for binomial, or conditional means for normal.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Additional arguments (unused).
#' @return An array with the same dimensions as \code{Y} (n1 x n2 x T)
#'   containing fitted values on the response scale.
#' @export
fitted.sir <- function(object, ...) {
  object$fitted.values
}

#' Extract Residuals from a SIR Model
#'
#' Returns residuals of the specified type. Response residuals are raw
#' (Y - fitted). Pearson residuals are standardized by the variance function.
#' Deviance residuals are signed square roots of the individual deviance
#' contributions and are most useful for diagnostic plots.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param type Character string specifying residual type: \code{"deviance"}
#'   (default), \code{"pearson"}, or \code{"response"}.
#' @param ... Additional arguments (unused).
#' @return An array with the same dimensions as \code{Y} containing the
#'   requested residuals. Contains NA where \code{Y} is missing.
#' @export
residuals.sir <- function(object, type = c("deviance", "pearson", "response"), ...) {
  type <- match.arg(type)
  if (is.null(object$residuals)) return(NULL)
  object$residuals[[type]]
}

#' Extract Log-Likelihood from a SIR Model
#'
#' Returns the log-likelihood at convergence as a \code{logLik} object, with
#' attributes for degrees of freedom and number of observations. This allows
#' \code{AIC()} and \code{BIC()} to work directly on the result.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Additional arguments (unused).
#' @return A \code{logLik} object with attributes \code{df} (number of
#'   estimated parameters) and \code{nobs} (number of observations).
#' @export
logLik.sir <- function(object, ...) {
  val <- object$ll
  attr(val, "df") <- length(object$summ$coef)
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}

#' Extract Number of Observations from a SIR Model
#'
#' Returns the number of non-missing dyad-time observations used in fitting.
#' For square networks, diagonal entries (self-loops) are excluded from this
#' count.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Additional arguments (unused).
#' @return Integer count of observations.
#' @export
nobs.sir <- function(object, ...) {
  object$nobs
}

#' Akaike Information Criterion for a SIR Model
#'
#' Computes AIC = -2 * log-likelihood + k * (number of parameters).
#' Use this to compare SIR models with different specifications (e.g.,
#' different numbers of influence covariates).
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Additional arguments for comparison with other models.
#' @param k Numeric penalty per parameter (default 2 for standard AIC).
#' @return Numeric AIC value. Lower is better.
#' @export
AIC.sir <- function(object, ..., k = 2) {
  ll <- logLik(object)
  -2 * as.numeric(ll) + k * attr(ll, "df")
}

#' Bayesian Information Criterion for a SIR Model
#'
#' Computes BIC = -2 * log-likelihood + log(nobs) * (number of parameters).
#' BIC penalizes model complexity more heavily than AIC for large samples.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Additional arguments for comparison with other models.
#' @return Numeric BIC value. Lower is better.
#' @export
BIC.sir <- function(object, ...) {
  ll <- logLik(object)
  nobs <- attr(ll, "nobs")
  -2 * as.numeric(ll) + log(nobs) * attr(ll, "df")
}

#' Summary of a Fitted SIR Model
#'
#' Produces a detailed summary of the fitted model including coefficient
#' estimates with standard errors and p-values, model fit statistics
#' (log-likelihood, AIC, BIC), convergence status, and summaries of the
#' estimated influence matrices A and B.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Additional arguments (unused).
#' @return An object of class \code{"summary.sir"} containing:
#'   \describe{
#'     \item{coefficients}{Data frame with columns \code{coef}, \code{se},
#'       \code{p.value}, and significance codes.}
#'     \item{loglik}{Log-likelihood at convergence.}
#'     \item{aic}{AIC value.}
#'     \item{bic}{BIC value.}
#'     \item{converged}{Logical convergence indicator.}
#'     \item{iterations}{Iteration count.}
#'     \item{A.summary}{List with mean, sd, and range of off-diagonal
#'       entries in the sender effects matrix.}
#'     \item{B.summary}{Same for the receiver effects matrix.}
#'   }
#' @seealso \code{\link{print.summary.sir}} for the printed output.
#' @export
summary.sir <- function(object, ...) {

  # create summary object
  ans <- list()
  ans$call <- object$call
  ans$family <- object$family
  ans$method <- object$method
  ans$symmetric <- isTRUE(object$symmetric)
  ans$fix_receiver <- isTRUE(object$fix_receiver)

  # model dimensions
  ans$m <- if (!is.null(object$m)) object$m else nrow(object$A)
  ans$n1 <- if (!is.null(object$n1)) object$n1 else nrow(object$A)
  ans$n2 <- if (!is.null(object$n2)) object$n2 else ncol(object$B)
  ans$bipartite <- isTRUE(object$bipartite)
  ans$T <- object$n_periods
  ans$p <- object$p
  ans$q <- object$q
  ans$nobs <- object$nobs

  # coefficients table with significance
  ans$coefficients <- object$summ
  if ("se" %in% colnames(ans$coefficients)) {
	z_scores <- ans$coefficients$coef / ans$coefficients$se
	ans$coefficients$p.value <- 2 * (1 - pnorm(abs(z_scores)))

	ans$coefficients$sig <- symnum(ans$coefficients$p.value,
									corr = FALSE, na = FALSE,
									cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
									symbols = c("***", "**", "*", ".", " "))
  }

  # model fit statistics
  ans$loglik <- object$ll
  ans$aic <- AIC(object)
  ans$bic <- BIC(object)
  ans$deviance <- if (!is.null(object$deviance)) object$deviance else -2 * object$ll
  ans$null.deviance <- object$null.deviance

  # convergence info
  ans$converged <- object$convergence
  ans$iterations <- object$iterations

  # influence matrices summary
  # handle dynamic W case where A/B are 3D arrays
  if (length(dim(object$A)) == 3) {
	# dynamic W: summarize across all time slices
	A_vals <- c()
	B_vals <- c()
	for (t in seq_len(dim(object$A)[3])) {
	  A_t <- object$A[,,t]
	  B_t <- object$B[,,t]
	  if (nrow(A_t) == ncol(A_t)) {
		A_vals <- c(A_vals, A_t[row(A_t) != col(A_t)])
	  } else {
		A_vals <- c(A_vals, c(A_t))
	  }
	  if (nrow(B_t) == ncol(B_t)) {
		B_vals <- c(B_vals, B_t[row(B_t) != col(B_t)])
	  } else {
		B_vals <- c(B_vals, c(B_t))
	  }
	}
	A_offdiag <- A_vals
	B_offdiag <- B_vals
  } else if (nrow(object$A) == ncol(object$A)) {
	A_offdiag <- object$A[row(object$A) != col(object$A)]
	B_offdiag <- if (nrow(object$B) == ncol(object$B)) {
	  object$B[row(object$B) != col(object$B)]
	} else {
	  c(object$B)
	}
  } else {
	A_offdiag <- c(object$A)
	B_offdiag <- c(object$B)
  }

  ans$A.summary <- list(
	mean = mean(A_offdiag),
	sd = sd(A_offdiag),
	range = range(A_offdiag)
  )

  ans$B.summary <- list(
	mean = mean(B_offdiag),
	sd = sd(B_offdiag),
	range = range(B_offdiag)
  )

  # sigma for normal family
  if (object$family == "normal" && !is.null(object$sigma2)) {
	ans$sigma <- sqrt(object$sigma2)
  }

  class(ans) <- "summary.sir"
  ans
}

#' Print a SIR Model Summary
#'
#' Displays the full model summary including network dimensions, family,
#' method, a coefficient table with significance stars (when standard errors
#' are available), model fit statistics, convergence status, and influence
#' matrix summaries.
#'
#' @param x A \code{summary.sir} object from \code{\link{summary.sir}}.
#' @param digits Number of significant digits to print. Default uses
#'   \code{getOption("digits") - 3}.
#' @param signif.stars Logical, whether to show significance stars beside
#'   p-values. Default uses \code{getOption("show.signif.stars")}.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the summary object.
#' @export
print.summary.sir <- function(x, digits = max(3L, getOption("digits") - 3L),
							   signif.stars = getOption("show.signif.stars"), ...) {

  cli::cli_h1("Social Influence Regression Model")

  # network and model info
  if (isTRUE(x$bipartite)) {
	net_type <- "bipartite"
	cli::cli_text("{.strong Network:} {.val {x$n1}} senders x {.val {x$n2}} receivers, {.val {x$n_periods}} time periods ({net_type})")
  } else {
	net_type <- if (isTRUE(x$symmetric)) "symmetric (undirected)" else "directed"
	cli::cli_text("{.strong Network:} {.val {x$m}} nodes, {.val {x$n_periods}} time periods ({net_type})")
  }
  cli::cli_text("{.strong Family:} {.val {x$family}} | {.strong Method:} {.val {x$method}}")
  if (isTRUE(x$fix_receiver) && !isTRUE(x$symmetric)) {
	cli::cli_text("{.strong Receiver:} fixed (B = I)")
  }
  if (!is.null(x$nobs)) {
	cli::cli_text("{.strong Observations:} {.val {x$nobs}}")
  }

  cli::cli_rule()

  # coefficients
  cli::cli_h2("Coefficients")

  if (nrow(x$coefficients) > 0) {
	if ("se" %in% colnames(x$coefficients) && "p.value" %in% colnames(x$coefficients)) {
	  coef_mat <- as.matrix(x$coefficients[, c("coef", "se", "t_se", "p.value")])
	  colnames(coef_mat) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

	  if (signif.stars && "sig" %in% colnames(x$coefficients)) {
		printCoefmat(coef_mat, digits = digits, signif.stars = TRUE,
					 P.values = TRUE, has.Pvalue = TRUE)
	  } else {
		print(round(coef_mat, digits))
	  }
	} else {
	  # no SEs available
	  coef_mat <- as.matrix(x$coefficients[, "coef", drop = FALSE])
	  colnames(coef_mat) <- "Estimate"
	  print(round(coef_mat, digits))
	  cli::cli_text("{.emph (Standard errors not computed. Use calc_se = TRUE or boot_sir() for inference.)}")
	}
  } else {
	cli::cli_text("{.emph No coefficients estimated.}")
  }

  cli::cli_rule()

  # model fit
  cli::cli_h2("Model Fit")

  fit_stats <- c(
	paste0("Log-Likelihood: ", sprintf("%.2f", x$loglik)),
	paste0("AIC: ", sprintf("%.2f", x$aic)),
	paste0("BIC: ", sprintf("%.2f", x$bic))
  )

  if (!is.null(x$sigma)) {
	fit_stats <- c(fit_stats,
				   paste0("Residual std. error: ", sprintf("%.4f", x$sigma)))
  }

  cli::cli_ul(fit_stats)

  # convergence
  if (x$converged) {
	cli::cli_alert_success("Converged in {.val {x$iterations}} iterations")
  } else {
	cli::cli_alert_warning("Did not converge")
  }

  cli::cli_rule()

  # influence matrices summary
  if (x$p > 0) {
	if (isTRUE(x$symmetric)) {
	  cli::cli_h2("Influence Matrix")
	  cli::cli_text("{.strong A matrix (influence weights):}")
	} else {
	  cli::cli_h2("Influence Matrices")
	  cli::cli_text("{.strong A matrix (sender effects):}")
	}

	cli::cli_ul(c(
	  paste0("Mean: ", sprintf("%.4f", x$A.summary$mean)),
	  paste0("SD: ", sprintf("%.4f", x$A.summary$sd)),
	  paste0("Range: [", sprintf("%.4f", x$A.summary$range[1]), ", ",
			 sprintf("%.4f", x$A.summary$range[2]), "]")
	))

	if (!isTRUE(x$fix_receiver)) {
	  cli::cli_text("{.strong B matrix (receiver effects):}")
	  cli::cli_ul(c(
		paste0("Mean: ", sprintf("%.4f", x$B.summary$mean)),
		paste0("SD: ", sprintf("%.4f", x$B.summary$sd)),
		paste0("Range: [", sprintf("%.4f", x$B.summary$range[1]), ", ",
			   sprintf("%.4f", x$B.summary$range[2]), "]")
	  ))
	}
  }

  invisible(x)
}


#' Print a Fitted SIR Model
#'
#' Displays a compact overview of the fitted model: network dimensions,
#' family, method, convergence status, log-likelihood, AIC, and coefficient
#' estimates. Use \code{summary()} for a more detailed report with p-values
#' and influence matrix summaries.
#'
#' @param x A fitted \code{sir} object from \code{\link{sir}}.
#' @param digits Number of digits to print. Default uses
#'   \code{getOption("digits") - 3}.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the sir object.
#' @export
print.sir <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cli::cli_text("\n")
  cli::cli_text("{.strong Social Influence Regression Model}")

  # network info
  if (isTRUE(x$bipartite)) {
	cli::cli_text("{.val {x$n1}} senders x {.val {x$n2}} receivers, {.val {x$n_periods}} time periods (bipartite)")
  } else {
	m <- if (!is.null(x$m)) x$m else nrow(x$A)
	net_type <- if (isTRUE(x$symmetric)) "symmetric" else "directed"
	cli::cli_text("{.val {m}} nodes, {.val {x$n_periods}} time periods ({net_type})")
  }

  # model config
  config_parts <- c(x$family, x$method)
  if (isTRUE(x$fix_receiver) && !isTRUE(x$symmetric)) {
	  config_parts <- c(config_parts, "fix_receiver")
  }
  cli::cli_text("Config: {.field {paste(config_parts, collapse = ' | ')}}")

  # status line
  n_obs <- if (!is.null(x$nobs)) x$nobs else NA
  if (x$convergence) {
	cli::cli_text("Status: {.field converged} | N = {.val {n_obs}} | Log-Lik: {.val {round(x$ll, 2)}} | AIC: {.val {round(AIC(x), 1)}}")
  } else {
	cli::cli_text("Status: {.emph not converged} | N = {.val {n_obs}} | Log-Lik: {.val {round(x$ll, 2)}}")
  }

  # coefficients as formatted table
  if (nrow(x$summ) > 0) {
	cli::cli_text("\nCoefficients:")
	if ("se" %in% colnames(x$summ)) {
	  coef_mat <- as.matrix(x$summ[, c("coef", "se")])
	  colnames(coef_mat) <- c("Estimate", "Std. Err")
	  print(round(coef_mat, digits))
	} else {
	  coef_mat <- as.matrix(x$summ[, "coef", drop = FALSE])
	  colnames(coef_mat) <- "Estimate"
	  print(round(coef_mat, digits))
	  cli::cli_text("{.emph (SEs not computed)}")
	}
  }

  cli::cli_text("\nUse {.code summary()} for detailed results")
  invisible(x)
}

#' Predictions from a Fitted SIR Model
#'
#' Generates predictions from a fitted SIR model for the training data or
#' for new data. Predictions can be on the link scale (linear predictor) or
#' the response scale (expected counts, probabilities, or means).
#'
#' For scenario (counterfactual) analysis, supply modified arrays in
#' \code{newdata}. For example, to see how the network would change if a
#' covariate increased by one unit, pass the modified Z array while keeping
#' W and X from the original fit.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param newdata Optional named list with components \code{W} (3D or 4D
#'   array), \code{X} (3D array), and/or \code{Z} (3D or 4D array) for
#'   counterfactual prediction. Dimensions must match the original fit. Any
#'   component not supplied is taken from the original fit. If NULL (default),
#'   returns predictions for the training data. Note: unlike many R predict
#'   methods, \code{newdata} is a list of arrays, not a data frame.
#' @param type Character string: \code{"link"} for linear predictor (eta)
#'   or \code{"response"} for expected values on the original scale.
#'   Default is \code{"response"}.
#' @param ... Additional arguments (unused).
#' @return An array (n1 x n2 x T) of predicted values on the requested scale.
#'
#' @examples
#' \dontrun{
#' model <- sir(Y, W, X, Z = Z, family = "poisson")
#'
#' # In-sample fitted values
#' pred <- predict(model)
#'
#' # Scenario: what if Z increases by 1 unit?
#' Z_shift <- Z + 1
#' pred_scenario <- predict(model, newdata = list(Z = Z_shift))
#'
#' # Compare mean predictions
#' mean(pred, na.rm = TRUE)
#' mean(pred_scenario, na.rm = TRUE)
#' }
#' @export
predict.sir <- function(object, newdata = NULL,
						type = c("response", "link"), ...) {
  type <- match.arg(type)

  fr <- isTRUE(object$fix_receiver)

  if (!is.null(newdata)) {
	# extract components from newdata, falling back to original fit
	W_new <- if (!is.null(newdata$W)) newdata$W else object$W
	X_new <- if (!is.null(newdata$X)) newdata$X else object$X
	Z_new <- if (!is.null(newdata$Z)) newdata$Z else object$Z

	# calculate linear predictor
	eta <- eta_tab(object$tab, W_new, X_new, Z_new, fix_receiver=fr)

	if (type == "response") {
	  if (object$family == "poisson") {
		return(exp(eta))
	  } else if (object$family == "binomial") {
		return(1 / (1 + exp(-eta)))
	  } else {
		return(eta)
	  }
	} else {
	  return(eta)
	}
  } else {
	# return predictions from training data
	if (type == "response") {
		return(object$fitted.values)
	} else {
		return(eta_tab(object$tab, object$W, object$X, object$Z, fix_receiver=fr))
	}
  }
}

#' Variance-Covariance Matrix for SIR Model Parameters
#'
#' Returns the variance-covariance matrix of the estimated parameters. Two
#' types are available: classical (inverse Hessian) and robust (sandwich
#' estimator). The robust version is more reliable when the model is
#' misspecified or when observations are not independent.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param type Character string: \code{"classical"} (default) for
#'   Hessian-based vcov, or \code{"robust"} for the sandwich estimator.
#' @param ... Additional arguments (unused).
#' @return A square matrix with rows and columns named by parameter.
#'   Returns NULL if standard errors were not computed (\code{calc_se = FALSE}).
#' @seealso \code{\link{confint.sir}} for confidence intervals.
#' @export
vcov.sir <- function(object, type = c("classical", "robust"), ...) {
  type <- match.arg(type)
  if (type == "robust") {
	  return(object$vcov_robust)
  }
  object$vcov
}

#' Confidence Intervals for SIR Model Parameters
#'
#' Computes confidence intervals using either Wald-based intervals (from the
#' Hessian standard errors) or bootstrap percentile intervals. Wald intervals
#' are the default and require that the model was fit with \code{calc_se = TRUE}.
#' Bootstrap intervals require a \code{\link{boot_sir}} result and tend to
#' be more reliable when the Hessian is ill-conditioned.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param parm Character vector of parameter names or integer indices to
#'   include. If NULL (default), returns intervals for all parameters.
#' @param level Confidence level between 0 and 1. Default is 0.95.
#' @param boot Optional \code{boot_sir} object from \code{\link{boot_sir}}.
#'   When provided, percentile intervals from the bootstrap distribution are
#'   used instead of Wald intervals.
#' @param ... Additional arguments (unused).
#' @return A matrix with one row per parameter and columns for the lower
#'   and upper bounds, labeled by percentage (e.g., \code{"2.5 \%"} and
#'   \code{"97.5 \%"}).
#' @seealso \code{\link{boot_sir}} for bootstrap inference,
#'   \code{\link{vcov.sir}} for the variance-covariance matrix.
#' @export
confint.sir <- function(object, parm = NULL, level = 0.95, boot = NULL, ...) {
  cf <- coef(object)
  pnames <- names(cf)
  a <- (1 - level) / 2
  pct <- paste0(format(100 * c(a, 1 - a), trim = TRUE, digits = 3), " %")

  if (!is.null(boot)) {
	# bootstrap percentile intervals
	valid <- apply(boot$coefs, 1, function(r) !any(is.na(r)))
	boot_valid <- boot$coefs[valid, , drop = FALSE]
	ci <- t(apply(boot_valid, 2, quantile, probs = c(a, 1 - a)))
	rownames(ci) <- pnames
	colnames(ci) <- pct
  } else {
	# wald intervals from vcov
	se <- object$summ$se
	if (is.null(se) || all(is.na(se))) {
	  cli::cli_abort("No standard errors available. Refit with {.code calc_se = TRUE} or provide bootstrap results via the {.arg boot} argument.")
	}
	z <- qnorm(1 - a)
	ci <- cbind(cf - z * se, cf + z * se)
	rownames(ci) <- pnames
	colnames(ci) <- pct
  }

  if (!is.null(parm)) {
	if (is.character(parm)) {
	  ci <- ci[parm, , drop = FALSE]
	} else {
	  ci <- ci[parm, , drop = FALSE]
	}
  }

  ci
}