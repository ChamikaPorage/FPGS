#' Estimate the ATE Using FPGS Methods
#'
#' Estimates the average treatment effect (ATE) using regression imputation
#' (RI), matching, inverse probability weighting (IPW), normalized inverse
#' probability weighting (NIPW), augmented inverse probability weighting
#' (AIPW), and stratification based on the Full Prognostic Score (FPGS).
#'
#' The FPGS may either be estimated using \code{fpgs()} or supplied as
#' externally estimated values of its two components.
#'
#' @param fit An optional object of class \code{"fpgs"} returned by
#'   \code{fpgs()}.
#' @param mu0 Optional numeric vector containing externally estimated values
#'   of the FPGS component under treatment level 0.
#' @param mu1 Optional numeric vector containing externally estimated values
#'   of the FPGS component under treatment level 1.
#' @param data Data frame containing the observed data. Required when
#'   \code{mu0} and \code{mu1} are supplied.
#' @param outcome Character string giving the name of the outcome variable.
#'   Required when \code{mu0} and \code{mu1} are supplied.
#' @param treatment Character string giving the name of the treatment
#'   variable. Required when \code{mu0} and \code{mu1} are supplied.
#' @param outcome_type Type of outcome: \code{"continuous"} or \code{"binary"}.
#'   Required when externally estimated FPGS components are supplied.
#' @param variance Character string specifying the variance estimation
#'   method. Options are \code{"none"}, \code{"analytic"}, or
#'   \code{"bootstrap"}. Default is \code{"none"}. The implemented variance
#'   estimators are available only when the FPGS is estimated using
#'   \code{fpgs()}.
#' @param B Number of bootstrap replications when
#'   \code{variance = "bootstrap"}. Default is 1000.
#' @param conf.level Confidence level for confidence intervals.
#'   Default is 0.95.
#' @param M Number of matches used by \code{fpgs_match()}. Default is 1.
#' @param n_bins Number of quantile-based bins used for each FPGS component
#'   in \code{fpgs_stratify()}. Default is 4.
#' @param min_obs Minimum recommended number of treated and untreated
#'   observations within each retained stratum. Default is 3.
#'
#' @return A data frame containing treatment-effect estimates. When analytic
#'   or bootstrap variance estimation is requested, estimated variances,
#'   standard errors, and confidence intervals are also returned.
#'
#' @details
#' Either \code{fit} or the pair \code{mu0} and \code{mu1} must be supplied,
#' but not both.
#'
#' When \code{fit} is supplied, it must be an object returned by
#' \code{fpgs()}. Both point estimation and the supported variance estimation
#' procedures are available.
#'
#' Alternatively, users may supply externally estimated FPGS components
#' through \code{mu0} and \code{mu1}. This allows the FPGS to be estimated
#' using outcome-modeling procedures not implemented directly in the package.
#' For externally estimated FPGS components, only point estimates are
#' available and \code{variance = "none"} must be used.
#'
#' If \code{variance = "none"}, only point estimates are returned.
#'
#' If \code{variance = "analytic"}, analytic variance estimators are used.
#' This option is available only for FPGS objects estimated using
#' \code{method = "linear_regression"} or
#' \code{method = "logistic_regression"} and returned by \code{fpgs()}.
#' Regression imputation and weighting estimators use sandwich variance
#' estimation, matching uses the variance estimate returned by
#' \code{Matching::Match()}, and stratification uses an analytic
#' within-stratum variance estimator.
#'
#' If \code{variance = "bootstrap"}, bootstrap variance estimates, standard
#' errors, and percentile confidence intervals are calculated by resampling
#' observations, refitting the FPGS, and recomputing the treatment-effect
#' estimator. Ordinary nonparametric bootstrap variance estimation is not
#' performed for the matching estimator.
#'
#' @export
fpgs_ate <- function(fit = NULL,
                     mu0 = NULL,
                     mu1 = NULL,
                     data = NULL,
                     outcome = NULL,
                     treatment = NULL,
                     outcome_type = NULL,
                     variance = c("none", "analytic", "bootstrap"),
                     B = 1000,
                     conf.level = 0.95,
                     M = 1,
                     n_bins = 4,
                     min_obs = 3) {

  # ------------------------------------------------------------
  # Match variance argument
  # ------------------------------------------------------------

  variance <- match.arg(variance)


  # ------------------------------------------------------------
  # Determine input type
  # ------------------------------------------------------------

  internal <- !is.null(fit)
  external <- !is.null(mu0) || !is.null(mu1)


  # Neither input route supplied
  if (!internal && !external) {
    stop(
      "Supply either an FPGS object using 'fit' or externally ",
      "estimated FPGS components using 'mu0' and 'mu1'."
    )
  }


  # Both input routes supplied
  if (internal && external) {
    stop(
      "Supply either 'fit' or externally estimated 'mu0' and 'mu1', ",
      "but not both."
    )
  }


  # ------------------------------------------------------------
  # Validate FPGS object returned by fpgs()
  # ------------------------------------------------------------

  if (internal) {

    if (!inherits(fit, "fpgs")) {
      stop("'fit' must be an object returned by fpgs().")
    }
  }


  # ------------------------------------------------------------
  # Validate externally estimated FPGS components
  # ------------------------------------------------------------

  if (external) {

    if (is.null(mu0) || is.null(mu1)) {
      stop(
        "Both 'mu0' and 'mu1' must be supplied when using ",
        "externally estimated FPGS components."
      )
    }

    if (!is.numeric(mu0) || !is.numeric(mu1)) {
      stop("'mu0' and 'mu1' must be numeric vectors.")
    }

    if (is.null(data)) {
      stop(
        "'data' must be supplied when using externally ",
        "estimated FPGS components."
      )
    }

    if (!is.data.frame(data)) {
      stop("'data' must be a data frame.")
    }

    if (is.null(outcome)) {
      stop(
        "'outcome' must be supplied when using externally ",
        "estimated FPGS components."
      )
    }

    if (is.null(treatment)) {
      stop(
        "'treatment' must be supplied when using externally ",
        "estimated FPGS components."
      )
    }

    if (is.null(outcome_type)) {
      stop(
        "'outcome_type' must be supplied when using externally ",
        "estimated FPGS components."
      )
    }

    outcome_type <- match.arg(
      outcome_type,
      choices = c("continuous", "binary")
    )

    if (!is.character(outcome) || length(outcome) != 1) {
      stop("'outcome' must be a single character string.")
    }

    if (!is.character(treatment) || length(treatment) != 1) {
      stop("'treatment' must be a single character string.")
    }

    if (!outcome %in% names(data)) {
      stop("'outcome' was not found in 'data'.")
    }

    if (!treatment %in% names(data)) {
      stop("'treatment' was not found in 'data'.")
    }

    if (length(mu0) != nrow(data) ||
        length(mu1) != nrow(data)) {
      stop(
        "'mu0' and 'mu1' must contain one value for each ",
        "observation in 'data'."
      )
    }

    if (anyNA(mu0) || anyNA(mu1)) {
      stop("'mu0' and 'mu1' must not contain missing values.")
    }


    # Treatment must contain 0 and 1
    Tr <- data[[treatment]]

    if (!all(Tr %in% c(0, 1))) {
      stop("'treatment' must be coded as 0 and 1.")
    }

    if (length(unique(Tr)) != 2) {
      stop("'treatment' must contain both treatment groups.")
    }


    # ----------------------------------------------------------
    # Variance estimation is not implemented for external FPGS
    # ----------------------------------------------------------

    if (variance != "none") {
      stop(
        "The implemented variance estimators are available only when ",
        "the FPGS is estimated using fpgs(). For externally estimated ",
        "FPGS components, use variance = 'none'."
      )
    }


    # ----------------------------------------------------------
    # Construct object for point estimation only
    # ----------------------------------------------------------

    fit <- list(
      data = data,
      outcome = outcome,
      treatment = treatment,
      mu0_hat = as.numeric(mu0),
      mu1_hat = as.numeric(mu1),
      method = "external",
      outcome_type = outcome_type
    )

    class(fit) <- "fpgs_external"
  }


  # ------------------------------------------------------------
  # Validate confidence level
  # ------------------------------------------------------------

  if (!is.numeric(conf.level) ||
      length(conf.level) != 1 ||
      conf.level <= 0 ||
      conf.level >= 1) {

    stop("'conf.level' must be a number between 0 and 1.")
  }


  # ------------------------------------------------------------
  # Validate analytic variance request
  # ------------------------------------------------------------

  if (internal &&
      variance == "analytic" &&
      !fit$method %in% c(
        "linear_regression",
        "logistic_regression"
      )) {

    stop(
      "Analytic variance estimation is only available for FPGS fits ",
      "estimated using method = 'linear_regression' or ",
      "'logistic_regression'. Use variance = 'bootstrap' for ",
      "method = 'random_forest'."
    )
  }


  # ------------------------------------------------------------
  # Treatment-effect point estimates
  # ------------------------------------------------------------

  ri <- fpgs_ri(fit)

  matching <- fpgs_match(
    fit,
    M = M
  )

  ipw <- fpgs_weight(
    fit,
    type = "ipw"
  )

  nipw <- fpgs_weight(
    fit,
    type = "nipw"
  )

  aipw <- fpgs_weight(
    fit,
    type = "aipw"
  )

  stratification <- fpgs_stratify(
    fit,
    n_bins = n_bins,
    min_obs = min_obs
  )


  estimator_names <- c(
    "RI",
    "Matching",
    "IPW",
    "NIPW",
    "AIPW",
    "Stratification"
  )


  estimates <- c(
    ri$estimate,
    matching$estimate,
    ipw$estimate,
    nipw$estimate,
    aipw$estimate,
    stratification$estimate
  )


  # ------------------------------------------------------------
  # Point estimates only
  # ------------------------------------------------------------

  if (variance == "none") {

    return(
      data.frame(
        estimator = estimator_names,
        estimate = estimates,
        row.names = NULL
      )
    )
  }


  # ------------------------------------------------------------
  # Analytic variance estimation
  # ------------------------------------------------------------

  if (variance == "analytic") {

    ri_var <- ri_variance(fit)

    ipw_var <- weight_variance(
      fit,
      type = "ipw"
    )

    nipw_var <- weight_variance(
      fit,
      type = "nipw"
    )

    aipw_var <- weight_variance(
      fit,
      type = "aipw"
    )


    variances <- c(
      ri_var$variance,
      matching$variance,
      ipw_var$variance,
      nipw_var$variance,
      aipw_var$variance,
      stratification$variance
    )


    standard_errors <- c(
      ri_var$se,
      matching$se,
      ipw_var$se,
      nipw_var$se,
      aipw_var$se,
      stratification$se
    )


    z_value <- stats::qnorm(
      1 - (1 - conf.level) / 2
    )


    ci_lower <- estimates -
      z_value * standard_errors

    ci_upper <- estimates +
      z_value * standard_errors


    return(
      data.frame(
        estimator = estimator_names,
        estimate = estimates,
        variance = variances,
        se = standard_errors,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        row.names = NULL
      )
    )
  }


  # ------------------------------------------------------------
  # Bootstrap variance estimation
  # ------------------------------------------------------------

  bootstrap_results <- list(

    RI = fpgs_bootstrap(
      fit,
      estimator = "ri",
      B = B,
      conf.level = conf.level
    ),

    IPW = fpgs_bootstrap(
      fit,
      estimator = "ipw",
      B = B,
      conf.level = conf.level
    ),

    NIPW = fpgs_bootstrap(
      fit,
      estimator = "nipw",
      B = B,
      conf.level = conf.level
    ),

    AIPW = fpgs_bootstrap(
      fit,
      estimator = "aipw",
      B = B,
      conf.level = conf.level
    ),

    Stratification = fpgs_bootstrap(
      fit,
      estimator = "stratification",
      B = B,
      conf.level = conf.level,
      n_bins = n_bins,
      min_obs = min_obs
    )
  )


  bootstrap_table <- do.call(
    rbind,
    lapply(
      names(bootstrap_results),
      function(name) {

        result <- bootstrap_results[[name]]

        data.frame(
          estimator = name,
          estimate = result$estimate,
          variance = result$variance,
          se = result$se,
          ci_lower = result$conf.int[1],
          ci_upper = result$conf.int[2]
        )
      }
    )
  )


  # ------------------------------------------------------------
  # Matching point estimate
  #
  # Ordinary nonparametric bootstrap variance estimation is not
  # performed for the matching estimator.
  # ------------------------------------------------------------

  matching_row <- data.frame(
    estimator = "Matching",
    estimate = matching$estimate,
    variance = NA_real_,
    se = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_
  )


  # ------------------------------------------------------------
  # Return estimators in a consistent order
  # ------------------------------------------------------------

  result <- rbind(
    bootstrap_table[
      bootstrap_table$estimator == "RI",
      ,
      drop = FALSE
    ],
    matching_row,
    bootstrap_table[
      bootstrap_table$estimator == "IPW",
      ,
      drop = FALSE
    ],
    bootstrap_table[
      bootstrap_table$estimator == "NIPW",
      ,
      drop = FALSE
    ],
    bootstrap_table[
      bootstrap_table$estimator == "AIPW",
      ,
      drop = FALSE
    ],
    bootstrap_table[
      bootstrap_table$estimator == "Stratification",
      ,
      drop = FALSE
    ]
  )

  rownames(result) <- NULL

  result
}
