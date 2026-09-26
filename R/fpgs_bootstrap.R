#' Bootstrap Variance Estimation for FPGS Estimators
#'
#' Computes bootstrap variance estimates, standard errors, and percentile
#' confidence intervals for FPGS treatment effect estimators.
#'
#' @param fit An object of class \code{"fpgs"} returned by \code{fpgs()}.
#' @param estimator Character string specifying the estimator. One of
#'   \code{"ri"}, \code{"ipw"}, \code{"nipw"}, \code{"aipw"},
#'   or \code{"stratification"}.
#' @param B Number of bootstrap replications. Default is 1000.
#' @param conf.level Confidence level for the percentile bootstrap interval.
#'   Default is 0.95.
#' @param seed Optional random seed.
#' @param verbose Logical. If \code{TRUE}, prints bootstrap progress.
#' @param ... Additional arguments passed to the treatment effect estimator,
#'   such as \code{n_bins} and \code{min_obs} for
#'   \code{fpgs_stratify()}.
#'
#' @return An object of class \code{"fpgs_bootstrap"} containing the point
#'   estimate, bootstrap variance, bootstrap standard error, percentile
#'   confidence interval, bootstrap estimates, and number of valid bootstrap
#'   replications.
#'
#' @details
#' In each bootstrap replication, observations are sampled with replacement.
#' The FPGS is then re-estimated using the same outcome type and estimation
#' method as in the original fit. The selected treatment effect estimator is
#' subsequently recomputed using the re-estimated FPGS.
#'
#' For random forest FPGS fits, the original number of folds, number of trees,
#' \code{mtry}, and minimum node size are reused when available.
#'
#' Bootstrap standard errors are calculated as the standard deviation of the
#' valid bootstrap estimates. Bootstrap variance is the squared bootstrap
#' standard error. Percentile confidence intervals are obtained from the
#' empirical bootstrap quantiles.
#'
#' Matching is not included in this bootstrap function.
#'
#' @export
fpgs_bootstrap <- function(
    fit,
    estimator = c(
      "ri",
      "ipw",
      "nipw",
      "aipw",
      "stratification"
    ),
    B = 1000,
    conf.level = 0.95,
    seed = NULL,
    verbose = TRUE,
    ...) {

  # ----------------------------------------------------------
  # Check estimator
  # ----------------------------------------------------------

  estimator <- match.arg(estimator)

  # ----------------------------------------------------------
  # Check FPGS object
  # ----------------------------------------------------------

  if (!inherits(fit, "fpgs")) {
    stop(
      "fit must be an object returned by fpgs().",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Check B
  # ----------------------------------------------------------

  if (!is.numeric(B) ||
      length(B) != 1 ||
      !is.finite(B) ||
      B < 2 ||
      B != as.integer(B)) {

    stop(
      "B must be an integer greater than or equal to 2.",
      call. = FALSE
    )
  }

  B <- as.integer(B)

  # ----------------------------------------------------------
  # Check confidence level
  # ----------------------------------------------------------

  if (!is.numeric(conf.level) ||
      length(conf.level) != 1 ||
      !is.finite(conf.level) ||
      conf.level <= 0 ||
      conf.level >= 1) {

    stop(
      "conf.level must be a number between 0 and 1.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Set random seed if supplied
  # ----------------------------------------------------------

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ----------------------------------------------------------
  # Original data
  # ----------------------------------------------------------

  data <- as.data.frame(fit$data)

  n <- nrow(data)

  if (n < 2) {
    stop(
      "The data must contain at least two observations.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Original point estimate
  # ----------------------------------------------------------

  estimate <- .fpgs_apply_estimator(
    fit = fit,
    estimator = estimator,
    ...
  )

  # ----------------------------------------------------------
  # Storage for bootstrap estimates
  # ----------------------------------------------------------

  bootstrap_estimates <- rep(
    NA_real_,
    B
  )

  # ----------------------------------------------------------
  # Bootstrap loop
  # ----------------------------------------------------------

  for (b in seq_len(B)) {

    if (
      verbose &&
      (
        b == 1 ||
        b %% 25 == 0 ||
        b == B
      )
    ) {

      message(
        "Bootstrap replication ",
        b,
        " of ",
        B
      )
    }

    # Sample observations with replacement

    bootstrap_index <- sample(
      seq_len(n),
      size = n,
      replace = TRUE
    )

    bootstrap_data <- data[
      bootstrap_index,
      ,
      drop = FALSE
    ]

    # Refit FPGS and recompute selected estimator

    bootstrap_estimates[b] <- tryCatch(

      {

        bootstrap_fit <- .fpgs_refit(
          fit = fit,
          data = bootstrap_data
        )

        .fpgs_apply_estimator(
          fit = bootstrap_fit,
          estimator = estimator,
          ...
        )

      },

      error = function(e) {

        if (verbose) {
          message(
            "Bootstrap replication ",
            b,
            " failed: ",
            conditionMessage(e)
          )
        }

        NA_real_
      }
    )
  }

  # ----------------------------------------------------------
  # Keep valid bootstrap estimates
  # ----------------------------------------------------------

  valid_estimates <- bootstrap_estimates[
    is.finite(bootstrap_estimates)
  ]

  B_valid <- length(valid_estimates)

  if (B_valid < 2) {

    stop(
      "Fewer than two valid bootstrap estimates were obtained.",
      call. = FALSE
    )
  }

  if (B_valid < 0.9 * B) {

    warning(
      "More than 10% of bootstrap replications failed. ",
      "Bootstrap variance estimates may be unreliable.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Bootstrap variance and standard error
  # ----------------------------------------------------------

  se <- stats::sd(
    valid_estimates
  )

  variance <- se^2

  # ----------------------------------------------------------
  # Percentile confidence interval
  # ----------------------------------------------------------

  alpha <- 1 - conf.level

  conf_int <- stats::quantile(
    valid_estimates,
    probs = c(
      alpha / 2,
      1 - alpha / 2
    ),
    na.rm = TRUE,
    names = FALSE
  )

  # ----------------------------------------------------------
  # Output
  # ----------------------------------------------------------

  out <- list(

    estimate = estimate,

    variance = variance,

    se = se,

    conf.int = conf_int,

    conf.level = conf.level,

    boot_estimates = valid_estimates,

    boot_estimates_all = bootstrap_estimates,

    B = B,

    B_valid = B_valid,

    estimator = estimator,

    outcome_type = fit$outcome_type,

    fpgs_method = fit$method,

    method = paste(
      "Bootstrap variance estimation for FPGS",
      estimator
    )
  )

  class(out) <- "fpgs_bootstrap"

  out
}


# ============================================================
# Internal helper:
# Refit FPGS in a bootstrap sample
# ============================================================

.fpgs_refit <- function(
    fit,
    data
) {

  # ----------------------------------------------------------
  # Linear-regression FPGS
  # ----------------------------------------------------------

  if (fit$method == "linear_regression") {

    return(
      fpgs(
        data = data,
        outcome = fit$outcome,
        treatment = fit$treatment,
        covariates = fit$covariates,
        outcome_type = "continuous",
        method = "linear_regression"
      )
    )
  }

  # ----------------------------------------------------------
  # Logistic-regression FPGS
  # ----------------------------------------------------------

  if (fit$method == "logistic_regression") {

    return(
      fpgs(
        data = data,
        outcome = fit$outcome,
        treatment = fit$treatment,
        covariates = fit$covariates,
        outcome_type = "binary",
        method = "logistic_regression"
      )
    )
  }

  # ----------------------------------------------------------
  # Random-forest FPGS
  # ----------------------------------------------------------

  if (fit$method == "random_forest") {

    folds <- if (!is.null(fit$folds)) {
      fit$folds
    } else {
      5
    }

    num.trees <- if (!is.null(fit$num.trees)) {
      fit$num.trees
    } else {
      500
    }

    mtry <- if (!is.null(fit$mtry)) {
      fit$mtry
    } else {
      NULL
    }

    min.node.size <- if (!is.null(fit$min.node.size)) {
      fit$min.node.size
    } else {
      5
    }

    return(
      fpgs(
        data = data,
        outcome = fit$outcome,
        treatment = fit$treatment,
        covariates = fit$covariates,
        outcome_type = fit$outcome_type,
        method = "random_forest",
        folds = folds,
        num.trees = num.trees,
        mtry = mtry,
        min.node.size = min.node.size
      )
    )
  }

  stop(
    "Unsupported FPGS estimation method.",
    call. = FALSE
  )
}


# ============================================================
# Internal helper:
# Apply one treatment-effect estimator
# ============================================================

.fpgs_apply_estimator <- function(
    fit,
    estimator,
    ...
) {

  # ----------------------------------------------------------
  # Regression imputation
  # ----------------------------------------------------------

  if (estimator == "ri") {

    return(
      fpgs_ri(
        fit
      )$estimate
    )
  }

  # ----------------------------------------------------------
  # Weighting
  # ----------------------------------------------------------

  if (
    estimator %in%
    c(
      "ipw",
      "nipw",
      "aipw"
    )
  ) {

    return(
      fpgs_weight(
        fit,
        type = estimator
      )$estimate
    )
  }

  # ----------------------------------------------------------
  # Stratification
  # ----------------------------------------------------------

  if (estimator == "stratification") {

    return(
      fpgs_stratify(
        fit,
        ...
      )$estimate
    )
  }

  stop(
    "Unknown estimator.",
    call. = FALSE
  )
}


# ============================================================
# Print method
# ============================================================

#' Print Method for FPGS Bootstrap Results
#'
#' Prints the bootstrap variance estimation results for an FPGS
#' treatment effect estimator.
#'
#' @param x An object of class \code{"fpgs_bootstrap"}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.fpgs_bootstrap <- function(
    x,
    ...
) {

  cat(
    "\n",
    x$method,
    "\n",
    sep = ""
  )

  cat(
    "----------------------------------\n"
  )

  cat(
    "Estimator:",
    x$estimator,
    "\n"
  )

  # User-friendly FPGS method label

  method_label <- switch(
    x$fpgs_method,
    linear_regression = "linear regression",
    logistic_regression = "logistic regression",
    random_forest = "random forest",
    x$fpgs_method
  )

  cat(
    "FPGS method:",
    method_label,
    "\n"
  )

  cat(
    "Outcome type:",
    x$outcome_type,
    "\n"
  )

  cat(
    "Estimate:",
    x$estimate,
    "\n"
  )

  cat(
    "Bootstrap variance:",
    x$variance,
    "\n"
  )

  cat(
    "Bootstrap SE:",
    x$se,
    "\n"
  )

  cat(
    paste0(
      100 * x$conf.level,
      "% CI: ["
    ),
    x$conf.int[1],
    ", ",
    x$conf.int[2],
    "]\n",
    sep = ""
  )

  cat(
    "Bootstrap replications:",
    x$B_valid,
    "valid out of",
    x$B,
    "\n"
  )

  invisible(x)
}


# ============================================================
# Summary method
# ============================================================

#' Summary Method for FPGS Bootstrap Results
#'
#' Summarizes bootstrap variance estimation results for an FPGS
#' treatment effect estimator.
#'
#' @param object An object of class \code{"fpgs_bootstrap"}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A data frame containing the bootstrap results.
#'
#' @export
summary.fpgs_bootstrap <- function(
    object,
    ...
) {

  data.frame(

    estimator = object$estimator,

    fpgs_method = object$fpgs_method,

    outcome_type = object$outcome_type,

    estimate = object$estimate,

    variance = object$variance,

    se = object$se,

    ci_lower = object$conf.int[1],

    ci_upper = object$conf.int[2],

    conf_level = object$conf.level,

    B = object$B,

    B_valid = object$B_valid
  )
}


# ============================================================
# Confidence interval method
# ============================================================

#' Confidence Interval for FPGS Bootstrap Results
#'
#' Computes a percentile bootstrap confidence interval for an FPGS
#' treatment effect estimate.
#'
#' @param object An object of class \code{"fpgs_bootstrap"}.
#' @param parm Currently ignored.
#' @param level Confidence level. Defaults to the level stored in
#'   \code{object}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A numeric vector containing the lower and upper confidence limits.
#'
#' @export
confint.fpgs_bootstrap <- function(
    object,
    parm = NULL,
    level = object$conf.level,
    ...
) {

  if (
    !is.numeric(level) ||
    length(level) != 1 ||
    !is.finite(level) ||
    level <= 0 ||
    level >= 1
  ) {

    stop(
      "level must be a number between 0 and 1.",
      call. = FALSE
    )
  }

  alpha <- 1 - level

  conf_int <- stats::quantile(
    object$boot_estimates,
    probs = c(
      alpha / 2,
      1 - alpha / 2
    ),
    na.rm = TRUE,
    names = FALSE
  )

  names(conf_int) <- c(
    "lower",
    "upper"
  )

  conf_int
}


# ============================================================
# Coefficient method
# ============================================================

#' Extract FPGS Bootstrap Estimate
#'
#' Extracts the estimated average treatment effect from an FPGS
#' bootstrap object.
#'
#' @param object An object of class \code{"fpgs_bootstrap"}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return The estimated average treatment effect.
#'
#' @export
coef.fpgs_bootstrap <- function(
    object,
    ...
) {

  object$estimate
}


# ============================================================
# Variance-covariance method
# ============================================================

#' Variance-Covariance Matrix for FPGS Bootstrap Estimate
#'
#' Returns the bootstrap variance estimate as a one-by-one
#' variance-covariance matrix.
#'
#' @param object An object of class \code{"fpgs_bootstrap"}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A one-by-one matrix containing the bootstrap variance estimate.
#'
#' @export
vcov.fpgs_bootstrap <- function(
    object,
    ...
) {

  variance_matrix <- matrix(
    object$variance,
    nrow = 1,
    ncol = 1
  )

  rownames(
    variance_matrix
  ) <- "ATE"

  colnames(
    variance_matrix
  ) <- "ATE"

  variance_matrix
}
