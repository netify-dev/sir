#' S3 Methods for SIR Objects
#' 
#' @description
#' Collection of S3 methods for objects of class "sir", providing standard model
#' interface functions including print, summary, plot, coef, fitted, residuals,
#' and model comparison tools.
#'
#' @importFrom cli cli_h1 cli_h2 cli_text cli_ul cli_alert_success cli_alert_warning cli_rule
#' @importFrom crayon blue green red yellow bold
#' @importFrom stats printCoefmat deviance quantile pnorm symnum AIC BIC logLik fitted residuals

#' Extract Model Coefficients
#' @param object A sir object
#' @param ... Additional arguments (unused)
#' @export
coef.sir <- function(object, ...) {
  object$summ$coef
}

#' Extract Fitted Values
#' @param object A sir object
#' @param ... Additional arguments (unused)
#' @export
fitted.sir <- function(object, ...) {
  # Store the data used in fitting for fitted value calculation
  if (!is.null(object$fitted.values)) {
    return(object$fitted.values)
  }
  
  # If not stored, return NULL with message
  cli::cli_alert_warning("Fitted values not stored in model object. Refit with `return.data = TRUE`")
  return(NULL)
}

#' Extract Model Residuals
#' @param object A sir object
#' @param type Type of residuals: "deviance", "pearson", or "response"
#' @param ... Additional arguments (unused)
#' @export
residuals.sir <- function(object, type = c("deviance", "pearson", "response"), ...) {
  type <- match.arg(type)
  
  if (!is.null(object$residuals)) {
    if (type == "deviance" && !is.null(object$residuals$deviance)) {
      return(object$residuals$deviance)
    } else if (type == "pearson" && !is.null(object$residuals$pearson)) {
      return(object$residuals$pearson)
    } else if (type == "response" && !is.null(object$residuals$response)) {
      return(object$residuals$response)
    }
  }
  
  cli::cli_alert_warning("Residuals not stored in model object. Refit with `return.data = TRUE`")
  return(NULL)
}

#' Extract Log-Likelihood
#' @param object A sir object
#' @param ... Additional arguments (unused)
#' @export
logLik.sir <- function(object, ...) {
  val <- object$ll
  attr(val, "df") <- length(object$summ$coef)
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}

#' Calculate AIC
#' @param object A sir object
#' @param ... Additional arguments for comparison with other models
#' @param k Penalty parameter (default 2 for AIC)
#' @export
AIC.sir <- function(object, ..., k = 2) {
  ll <- logLik(object)
  -2 * as.numeric(ll) + k * attr(ll, "df")
}

#' Calculate BIC
#' @param object A sir object
#' @param ... Additional arguments for comparison with other models
#' @export
BIC.sir <- function(object, ...) {
  ll <- logLik(object)
  nobs <- attr(ll, "nobs")
  -2 * as.numeric(ll) + log(nobs) * attr(ll, "df")
}

#' Summary Method for SIR Objects
#' 
#' @param object A sir object
#' @param ... Additional arguments
#' @return An object of class "summary.sir"
#' @export
summary.sir <- function(object, ...) {
  
  # Create summary object
  ans <- list()
  ans$call <- object$call
  ans$family <- object$family
  ans$method <- object$method
  
  # Model dimensions
  ans$m <- nrow(object$A)  # Number of nodes
  ans$T <- object$T         # Number of time periods (if stored)
  ans$p <- object$p         # Number of influence covariates
  ans$q <- object$q         # Number of exogenous covariates
  
  # Coefficients table with significance
  ans$coefficients <- object$summ
  if ("se" %in% colnames(ans$coefficients)) {
    # Add p-values (two-tailed test)
    z_scores <- ans$coefficients$coef / ans$coefficients$se
    ans$coefficients$p.value <- 2 * (1 - pnorm(abs(z_scores)))
    
    # Add significance stars
    ans$coefficients$sig <- symnum(ans$coefficients$p.value, 
                                    corr = FALSE, na = FALSE,
                                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                    symbols = c("***", "**", "*", ".", " "))
  }
  
  # Model fit statistics
  ans$loglik <- object$ll
  ans$aic <- AIC(object)
  ans$bic <- BIC(object)
  ans$deviance <- if (!is.null(object$deviance)) object$deviance else -2 * object$ll
  ans$null.deviance <- object$null.deviance
  
  # Convergence info
  ans$converged <- object$convergence
  ans$iterations <- object$iterations
  
  # Influence matrices summary
  ans$A.summary <- list(
    mean = mean(object$A[row(object$A) != col(object$A)]),
    sd = sd(object$A[row(object$A) != col(object$A)]),
    range = range(object$A[row(object$A) != col(object$A)])
  )
  
  ans$B.summary <- list(
    mean = mean(object$B[row(object$B) != col(object$B)]),
    sd = sd(object$B[row(object$B) != col(object$B)]),
    range = range(object$B[row(object$B) != col(object$B)])
  )
  
  # Pseudo R-squared for non-normal families
  if (object$family != "normal" && !is.null(ans$null.deviance)) {
    ans$pseudo.r.squared <- 1 - (ans$deviance / ans$null.deviance)
  } else if (object$family == "normal" && !is.null(object$sigma2)) {
    ans$sigma <- sqrt(object$sigma2)
    # Could add R-squared if we have TSS
  }
  
  class(ans) <- "summary.sir"
  ans
}

#' Print Summary of SIR Model
#' @param x A summary.sir object
#' @param digits Number of significant digits to print
#' @param signif.stars Logical, whether to show significance stars
#' @param ... Additional arguments (unused)
#' @export
print.summary.sir <- function(x, digits = max(3L, getOption("digits") - 3L), 
                               signif.stars = getOption("show.signif.stars"), ...) {
  
  cli::cli_h1("Social Influence Regression Model Summary")
  
  # Model info
  cli::cli_text("{.strong Family:} {.val {x$family}}")
  cli::cli_text("{.strong Method:} {.val {x$method}}")
  if (!is.null(x$m)) {
    cli::cli_text("{.strong Network size:} {.val {x$m}} nodes")
  }
  if (!is.null(x$T)) {
    cli::cli_text("{.strong Time periods:} {.val {x$T}}")
  }
  
  cli::cli_rule()
  
  # Coefficients
  cli::cli_h2("Coefficients")
  
  if (nrow(x$coefficients) > 0) {
    coef_mat <- as.matrix(x$coefficients[, c("coef", "se", "t_se", "p.value")])
    colnames(coef_mat) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    
    if (signif.stars && "sig" %in% colnames(x$coefficients)) {
      # Print with significance stars
      printCoefmat(coef_mat, digits = digits, signif.stars = TRUE, 
                   P.values = TRUE, has.Pvalue = TRUE)
    } else {
      print(round(coef_mat, digits))
    }
  } else {
    cli::cli_text("{.emph No coefficients (intercept-only model)}")
  }
  
  cli::cli_rule()
  
  # Model fit
  cli::cli_h2("Model Fit")
  
  fit_stats <- c(
    paste0("Log-Likelihood: ", sprintf("%.2f", x$loglik)),
    paste0("AIC: ", sprintf("%.2f", x$aic)),
    paste0("BIC: ", sprintf("%.2f", x$bic))
  )
  
  if (!is.null(x$pseudo.r.squared)) {
    fit_stats <- c(fit_stats, 
                   paste0("Pseudo R-squared: ", sprintf("%.4f", x$pseudo.r.squared)))
  }
  
  if (!is.null(x$sigma)) {
    fit_stats <- c(fit_stats,
                   paste0("Residual std error: ", sprintf("%.4f", x$sigma)))
  }
  
  cli::cli_ul(fit_stats)
  
  # Convergence
  if (x$converged) {
    cli::cli_alert_success("Converged in {.val {x$iterations}} iterations")
  } else {
    cli::cli_alert_warning("Did not converge")
  }
  
  cli::cli_rule()
  
  # Influence matrices summary
  cli::cli_h2("Influence Matrices")
  
  cli::cli_text("{.strong A matrix (sender effects):}")
  cli::cli_ul(c(
    paste0("Mean: ", sprintf("%.4f", x$A.summary$mean)),
    paste0("SD: ", sprintf("%.4f", x$A.summary$sd)),
    paste0("Range: [", sprintf("%.4f", x$A.summary$range[1]), ", ",
           sprintf("%.4f", x$A.summary$range[2]), "]")
  ))
  
  cli::cli_text("{.strong B matrix (receiver effects):}")
  cli::cli_ul(c(
    paste0("Mean: ", sprintf("%.4f", x$B.summary$mean)),
    paste0("SD: ", sprintf("%.4f", x$B.summary$sd)),
    paste0("Range: [", sprintf("%.4f", x$B.summary$range[1]), ", ",
           sprintf("%.4f", x$B.summary$range[2]), "]")
  ))
  
  invisible(x)
}


#' Simplified Print Method for SIR Objects
#' 
#' @param x A sir object
#' @param digits Number of digits to print
#' @param ... Additional arguments
#' @export
print.sir <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cli::cli_text("\n")
  cli::cli_text("{.strong Social Influence Regression Model}")
  cli::cli_text("Family: {.field {x$family}} | Method: {.field {x$method}}")
  
  # Basic info
  if (x$convergence) {
    cli::cli_text("Status: {crayon::green('Converged')} | Log-Lik: {.val {round(x$ll, 2)}}")
  } else {
    cli::cli_text("Status: {crayon::yellow('Not converged')} | Log-Lik: {.val {round(x$ll, 2)}}")
  }
  
  # Coefficients
  if (nrow(x$summ) > 0) {
    cli::cli_text("\nCoefficients:")
    print(round(x$summ$coef, digits), ...)
  } else {
    cli::cli_text("\n{.emph No coefficients (intercept-only model)}")
  }
  
  cli::cli_text("\nUse {.code summary()} for detailed results")
  invisible(x)
}

#' Predict Method for SIR Objects
#' 
#' @param object A sir object
#' @param newdata New data for prediction
#' @param type Type of prediction
#' @param ... Additional arguments
#' @export
predict.sir <- function(object, newdata = NULL, 
                        type = c("link", "response"), ...) {
  type <- match.arg(type)
  
  if (!is.null(newdata)) {
    # Extract components from newdata
    Y_new <- newdata$Y
    W_new <- if (!is.null(newdata$W)) newdata$W else object$W
    X_new <- if (!is.null(newdata$X)) newdata$X else newdata$Y  # Use Y as X if not provided
    Z_new <- newdata$Z
    
    # Calculate linear predictor
    eta <- eta_tab(object$tab, W_new, X_new, Z_new)
    
    if (type == "response") {
      # Transform based on family
      if (object$family == "poisson") {
        return(exp(eta))
      } else if (object$family == "binomial") {
        return(1 / (1 + exp(-eta)))
      } else {
        return(eta)  # Normal family
      }
    } else {
      return(eta)  # Link scale
    }
  } else {
    cli::cli_alert_warning("Prediction requires newdata")
    return(NULL)
  }
}

#' Extract Variance-Covariance Matrix
#' @param object A sir object
#' @param ... Additional arguments (unused)
#' @export
vcov.sir <- function(object, ...) {
  if (!is.null(object$vcov)) {
    return(object$vcov)
  }
  
  # If not stored, try to compute from Hessian
  if (!is.null(object$hessian)) {
    vcov_mat <- solve(object$hessian)
    return(vcov_mat)
  }
  
  cli::cli_alert_warning("Variance-covariance matrix not available")
  return(NULL)
}