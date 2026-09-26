#' FPGS Weighting Estimator
#'
#' Estimates the average treatment effect (ATE) using weighting methods based
#' on propensity scores estimated from the two-dimensional Full Prognostic
#' Score (FPGS).
#'
#' @param fit An object of class \code{"fpgs"} returned by \code{fpgs()},
#'   or an object of class \code{"fpgs_external"} containing externally
#'   estimated FPGS components.
#' @param type Character string specifying the weighting estimator:
#'   \code{"ipw"}, \code{"nipw"}, or \code{"aipw"}.
#'
#' @return An object of class \code{"fpgs_weight"} containing the ATE estimate,
#'   treatment-specific means, estimated propensity scores, and fitted models.
#'
#' @details
#' The propensity score is estimated using logistic regression with the two
#' estimated FPGS components as predictors. The FPGS components may be
#' estimated using \code{fpgs()} or supplied externally.
#'
#' For AIPW, additional treatment-specific outcome models are fitted using
#' the estimated FPGS components as predictors. Logistic regression is used
#' for binary outcomes and linear regression is used for continuous outcomes.
#'
#' This function computes point estimates for both internally and externally
#' estimated FPGS components. The variance estimators currently implemented
#' in \code{FPGS} are not available for externally estimated FPGS components.
#'
#' @export
fpgs_weight <- function(fit, type = c("ipw", "nipw", "aipw")) {

  type <- match.arg(type)

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
  # Extract data
  # ----------------------------------------------------------

  data <- fit$data
  Y <- data[[fit$outcome]]
  Tr <- data[[fit$treatment]]

  # ----------------------------------------------------------
  # Check treatment
  # ----------------------------------------------------------

  if (!all(Tr %in% c(0, 1))) {
    stop(
      "Treatment must be coded as 0 and 1.",
      call. = FALSE
    )
  }

  if (!all(c(0, 1) %in% Tr)) {
    stop(
      "Both treatment groups must be present.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Check FPGS components
  # ----------------------------------------------------------

  if (anyNA(fit$mu0_hat) || anyNA(fit$mu1_hat)) {
    stop(
      "The estimated FPGS components must not contain missing values.",
      call. = FALSE
    )
  }

  fpgs_data <- data.frame(
    Y = Y,
    Tr = Tr,
    mu0_hat = fit$mu0_hat,
    mu1_hat = fit$mu1_hat
  )

  # ----------------------------------------------------------
  # FPGS-based propensity score model
  # ----------------------------------------------------------

  ps_fit <- stats::glm(
    Tr ~ mu0_hat + mu1_hat,
    data = fpgs_data,
    family = stats::binomial()
  )

  e_hat <- stats::predict(
    ps_fit,
    type = "response"
  )

  # Avoid propensity scores equal to exactly 0 or 1
  e_hat <- pmin(
    pmax(e_hat, 1e-6),
    1 - 1e-6
  )

  # ----------------------------------------------------------
  # Initialize outcome-model objects
  # ----------------------------------------------------------

  outcome_model1 <- NULL
  outcome_model0 <- NULL
  outcome_mean1 <- NULL
  outcome_mean0 <- NULL

  # ----------------------------------------------------------
  # IPW
  # ----------------------------------------------------------

  if (type == "ipw") {

    mu1 <- mean(
      Tr * Y / e_hat
    )

    mu0 <- mean(
      (1 - Tr) * Y / (1 - e_hat)
    )
  }

  # ----------------------------------------------------------
  # Normalized IPW
  # ----------------------------------------------------------

  if (type == "nipw") {

    mu1 <- sum(
      Tr * Y / e_hat
    ) /
      sum(
        Tr / e_hat
      )

    mu0 <- sum(
      (1 - Tr) * Y / (1 - e_hat)
    ) /
      sum(
        (1 - Tr) / (1 - e_hat)
      )
  }

  # ----------------------------------------------------------
  # AIPW
  # ----------------------------------------------------------

  if (type == "aipw") {

    # --------------------------------------------------------
    # Check outcome type
    # --------------------------------------------------------

    if (is.null(fit$outcome_type) ||
        !fit$outcome_type %in% c("continuous", "binary")) {

      stop(
        "outcome_type must be either 'continuous' or 'binary' ",
        "for the AIPW estimator.",
        call. = FALSE
      )
    }

    # --------------------------------------------------------
    # Binary outcome
    # --------------------------------------------------------

    if (fit$outcome_type == "binary") {

      outcome_model1 <- stats::glm(
        Y ~ mu0_hat + mu1_hat,
        data = fpgs_data,
        subset = Tr == 1,
        family = stats::binomial()
      )

      outcome_model0 <- stats::glm(
        Y ~ mu0_hat + mu1_hat,
        data = fpgs_data,
        subset = Tr == 0,
        family = stats::binomial()
      )

      outcome_mean1 <- stats::predict(
        outcome_model1,
        newdata = fpgs_data,
        type = "response"
      )

      outcome_mean0 <- stats::predict(
        outcome_model0,
        newdata = fpgs_data,
        type = "response"
      )
    }

    # --------------------------------------------------------
    # Continuous outcome
    # --------------------------------------------------------

    if (fit$outcome_type == "continuous") {

      outcome_model1 <- stats::lm(
        Y ~ mu0_hat + mu1_hat,
        data = fpgs_data,
        subset = Tr == 1
      )

      outcome_model0 <- stats::lm(
        Y ~ mu0_hat + mu1_hat,
        data = fpgs_data,
        subset = Tr == 0
      )

      outcome_mean1 <- stats::predict(
        outcome_model1,
        newdata = fpgs_data
      )

      outcome_mean0 <- stats::predict(
        outcome_model0,
        newdata = fpgs_data
      )
    }

    # --------------------------------------------------------
    # AIPW treatment-specific means
    # --------------------------------------------------------

    mu1 <- mean(
      outcome_mean1 +
        Tr * (Y - outcome_mean1) / e_hat
    )

    mu0 <- mean(
      outcome_mean0 +
        (1 - Tr) * (Y - outcome_mean0) / (1 - e_hat)
    )
  }

  # ----------------------------------------------------------
  # Average treatment effect
  # ----------------------------------------------------------

  ate <- mu1 - mu0

  # ----------------------------------------------------------
  # Output
  # ----------------------------------------------------------

  out <- list(
    estimate = ate,
    mu1 = mu1,
    mu0 = mu0,
    propensity_score = e_hat,
    ps_model = ps_fit,
    outcome_model1 = outcome_model1,
    outcome_model0 = outcome_model0,
    outcome_mean1 = outcome_mean1,
    outcome_mean0 = outcome_mean0,
    type = type,
    method = paste("FPGS", type, "weighting")
  )

  class(out) <- "fpgs_weight"

  out
}
