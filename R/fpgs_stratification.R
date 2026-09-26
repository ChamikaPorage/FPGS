#' FPGS Stratification Estimator
#'
#' Estimates the average treatment effect (ATE) by stratifying on the
#' two-dimensional estimated Full Prognostic Score (FPGS).
#'
#' Strata containing only treated or only control individuals are excluded,
#' and the remaining stratum weights are normalized to sum to one.
#'
#' @param fit An object of class \code{"fpgs"} returned by \code{fpgs()},
#'   or an object of class \code{"fpgs_external"} containing externally
#'   estimated FPGS components.
#' @param n_bins Number of quantile-based bins used for each FPGS component.
#'   Must be an integer greater than or equal to 2. Default is 4.
#' @param min_obs Minimum recommended number of treated and control
#'   observations within each retained stratum. Must be a positive integer.
#'   Default is 3.
#'
#' @return An object of class \code{"fpgs_stratify"} containing the ATE
#'   estimate, analytic variance and standard error when available, stratum
#'   assignments, stratum-specific effects, and information about excluded
#'   strata.
#'
#' @details
#' Each component of the estimated FPGS is divided into quantile-based bins,
#' and the final strata are formed by cross-classifying the two sets of bins.
#'
#' Within each retained stratum, the treatment effect is estimated as the
#' difference between the mean outcome among treated and control individuals.
#' Strata containing only treated or only control individuals are excluded.
#'
#' The overall treatment effect is calculated as a weighted average of the
#' retained stratum-specific treatment effects. The weights are proportional
#' to the stratum sizes and are normalized over the retained strata.
#'
#' For FPGS objects returned by \code{fpgs()}, the analytic variance is
#' calculated from the within-stratum treated and control outcome variances,
#' treating the estimated strata as fixed. At least two treated and two
#' control observations are required within each retained stratum to calculate
#' the analytic variance. If this condition is not satisfied, the treatment
#' effect is still returned, but the variance and standard error are set to
#' \code{NA}.
#'
#' For externally estimated FPGS components, the stratification estimator is
#' available for point estimation, but the variance estimators currently
#' implemented in \code{FPGS} are not available.
#'
#' Bootstrap variance estimation for FPGS objects returned by \code{fpgs()}
#' may be used to account for uncertainty associated with estimation of the
#' FPGS and construction of the strata.
#'
#' The \code{min_obs} argument is used as a diagnostic only. Strata are not
#' excluded solely because they contain fewer than \code{min_obs} treated or
#' control observations.
#'
#' @export
fpgs_stratify <- function(fit, n_bins = 4, min_obs = 3) {

  # ----------------------------------------------------------
  # Check input object
  # ----------------------------------------------------------

  if (!inherits(fit, c("fpgs", "fpgs_external"))) {
    stop(
      "fit must be an object of class 'fpgs' or 'fpgs_external'.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Check n_bins
  # ----------------------------------------------------------

  if (!is.numeric(n_bins) ||
      length(n_bins) != 1 ||
      !is.finite(n_bins) ||
      n_bins < 2 ||
      n_bins != as.integer(n_bins)) {

    stop(
      "n_bins must be an integer greater than or equal to 2.",
      call. = FALSE
    )
  }

  n_bins <- as.integer(n_bins)

  # ----------------------------------------------------------
  # Check min_obs
  # ----------------------------------------------------------

  if (!is.numeric(min_obs) ||
      length(min_obs) != 1 ||
      !is.finite(min_obs) ||
      min_obs < 1 ||
      min_obs != as.integer(min_obs)) {

    stop(
      "min_obs must be a positive integer.",
      call. = FALSE
    )
  }

  min_obs <- as.integer(min_obs)

  # ----------------------------------------------------------
  # Extract data
  # ----------------------------------------------------------

  data <- fit$data
  Y <- data[[fit$outcome]]
  Tr <- data[[fit$treatment]]

  # ----------------------------------------------------------
  # Check FPGS components
  # ----------------------------------------------------------

  if (anyNA(fit$mu0_hat) || anyNA(fit$mu1_hat)) {
    stop(
      "The estimated FPGS components must not contain missing values.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Divide each FPGS component into quantile-based bins
  # ----------------------------------------------------------

  breaks_mu0 <- unique(
    stats::quantile(
      fit$mu0_hat,
      probs = seq(0, 1, length.out = n_bins + 1),
      na.rm = TRUE
    )
  )

  breaks_mu1 <- unique(
    stats::quantile(
      fit$mu1_hat,
      probs = seq(0, 1, length.out = n_bins + 1),
      na.rm = TRUE
    )
  )

  # Need at least two unique break points
  if (length(breaks_mu0) < 2 ||
      length(breaks_mu1) < 2) {

    stop(
      "The estimated FPGS components do not contain enough unique ",
      "values to construct strata.",
      call. = FALSE
    )
  }

  bin_mu0 <- cut(
    fit$mu0_hat,
    breaks = breaks_mu0,
    include.lowest = TRUE,
    labels = FALSE
  )

  bin_mu1 <- cut(
    fit$mu1_hat,
    breaks = breaks_mu1,
    include.lowest = TRUE,
    labels = FALSE
  )

  # ----------------------------------------------------------
  # Cross-classify the two FPGS components
  # ----------------------------------------------------------

  strata <- interaction(
    bin_mu0,
    bin_mu1,
    drop = TRUE
  )

  stratification_data <- data.frame(
    Y = Y,
    Tr = Tr,
    strata = strata
  )

  # ----------------------------------------------------------
  # Treatment effect within each stratum
  # ----------------------------------------------------------

  stratum_effects <- lapply(
    split(
      stratification_data,
      stratification_data$strata
    ),
    function(stratum_data) {

      n <- nrow(stratum_data)

      n_treated <- sum(
        stratum_data$Tr == 1
      )

      n_control <- sum(
        stratum_data$Tr == 0
      )

      # Treatment effect cannot be estimated without both groups
      if (n_treated == 0 || n_control == 0) {

        return(
          data.frame(
            n = n,
            n_treated = n_treated,
            n_control = n_control,
            mean_treated = NA_real_,
            mean_control = NA_real_,
            te = NA_real_
          )
        )
      }

      mean_treated <- mean(
        stratum_data$Y[
          stratum_data$Tr == 1
        ]
      )

      mean_control <- mean(
        stratum_data$Y[
          stratum_data$Tr == 0
        ]
      )

      data.frame(
        n = n,
        n_treated = n_treated,
        n_control = n_control,
        mean_treated = mean_treated,
        mean_control = mean_control,
        te = mean_treated - mean_control
      )
    }
  )

  stratum_effects <- do.call(
    rbind,
    stratum_effects
  )

  # ----------------------------------------------------------
  # Exclude strata without treatment overlap
  # ----------------------------------------------------------

  valid_strata <- !is.na(
    stratum_effects$te
  )

  n_dropped_strata <- sum(
    !valid_strata
  )

  n_dropped_obs <- sum(
    stratum_effects$n[
      !valid_strata
    ]
  )

  retained <- stratum_effects[
    valid_strata,
    ,
    drop = FALSE
  ]

  if (nrow(retained) == 0) {
    stop(
      "No strata contain both treated and control individuals.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Renormalize weights over retained strata
  # ----------------------------------------------------------

  n_retained_obs <- sum(
    retained$n
  )

  retained$weight <-
    retained$n / n_retained_obs

  retained$weighted_te <-
    retained$weight * retained$te

  # ----------------------------------------------------------
  # Flag small cells
  # ----------------------------------------------------------

  small_cells <- with(
    retained,
    n_treated < min_obs |
      n_control < min_obs
  )

  if (any(small_cells)) {
    warning(
      "Some retained strata contain fewer than ",
      min_obs,
      " treated or control observations.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Report strata without treatment overlap
  # ----------------------------------------------------------

  if (n_dropped_strata > 0) {
    warning(
      n_dropped_strata,
      " strata containing ",
      n_dropped_obs,
      " observations were excluded because they lacked treatment overlap.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Overall treatment effect
  # ----------------------------------------------------------

  ate <- sum(
    retained$weighted_te
  )

  # ----------------------------------------------------------
  # Analytic variance
  #
  # Only provided for FPGS estimated by fpgs()
  # ----------------------------------------------------------

  variance <- NA_real_
  se <- NA_real_

  if (inherits(fit, "fpgs")) {

    variance_components <- numeric(
      nrow(retained)
    )

    variance_available <- TRUE

    retained_names <- rownames(
      retained
    )

    for (j in seq_len(nrow(retained))) {

      stratum_name <- retained_names[j]

      in_stratum <-
        as.character(strata) == stratum_name

      Y_treated <- Y[
        in_stratum & Tr == 1
      ]

      Y_control <- Y[
        in_stratum & Tr == 0
      ]

      n_treated <- length(
        Y_treated
      )

      n_control <- length(
        Y_control
      )

      # Sample variance requires at least two observations
      # in each treatment group
      if (n_treated < 2 ||
          n_control < 2) {

        variance_available <- FALSE
        break
      }

      variance_treated <- stats::var(
        Y_treated
      )

      variance_control <- stats::var(
        Y_control
      )

      variance_components[j] <-
        retained$weight[j]^2 *
        (
          variance_treated / n_treated +
            variance_control / n_control
        )
    }

    if (variance_available) {

      variance <- sum(
        variance_components
      )

      se <- sqrt(
        variance
      )

    } else {

      warning(
        "Analytic variance could not be calculated because at least one ",
        "retained stratum contains fewer than two treated or control ",
        "observations.",
        call. = FALSE
      )
    }
  }

  # ----------------------------------------------------------
  # Output
  # ----------------------------------------------------------

  out <- list(
    estimate = ate,
    variance = variance,
    se = se,
    n_bins = n_bins,
    min_obs = min_obs,
    strata = strata,
    stratum_effects = stratum_effects,
    retained_stratum_effects = retained,
    n_retained_strata = nrow(retained),
    n_dropped_strata = n_dropped_strata,
    n_dropped_obs = n_dropped_obs,
    n_retained_obs = n_retained_obs,
    method = "FPGS stratification"
  )

  class(out) <- "fpgs_stratify"

  out
}
