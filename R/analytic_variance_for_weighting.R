#' Analytic Variance for FPGS Weighting Estimators
#'
#' Computes sandwich variance estimates and standard errors for the
#' FPGS inverse probability weighting (IPW), normalized inverse probability
#' weighting (NIPW), and augmented inverse probability weighting (AIPW)
#' estimators using stacked estimating equations.
#'
#' Both continuous and binary outcomes are supported. For continuous
#' outcomes, linear outcome models are used. For binary outcomes,
#' logistic outcome models are used.
#'
#' @param fit An object of class \code{"fpgs"} returned by \code{fpgs()}
#'   with \code{method = "linear_regression"} or
#'   \code{method = "logistic_regression"}.
#' @param type Character string specifying the weighting estimator.
#'   Must be one of \code{"ipw"}, \code{"nipw"}, or \code{"aipw"}.
#'
#' @return A list containing:
#' \describe{
#'   \item{estimate}{Estimated average treatment effect.}
#'   \item{variance}{Estimated sandwich variance of the ATE estimator.}
#'   \item{se}{Estimated standard error of the ATE estimator.}
#' }
#'
#' @details
#' The variance estimator is based on stacked estimating equations that
#' account for estimation of the two treatment-specific regression models
#' used to construct the FPGS and the propensity score model fitted using
#' the estimated two-dimensional FPGS.
#'
#' For AIPW, additional treatment-specific outcome models are fitted using
#' the two estimated FPGS components as predictors. The corresponding
#' estimating equations are included in the stacked system.
#'
#' This analytic variance estimator is available when the FPGS is estimated
#' using linear regression or logistic regression. Bootstrap variance
#' estimation may be used when the FPGS is estimated using random forest.
#'
#' @keywords internal
weight_variance <- function(fit,
                            type = c("ipw", "nipw", "aipw")) {

  # ------------------------------------------------------------
  # Check input
  # ------------------------------------------------------------

  if (!inherits(fit, "fpgs")) {
    stop(
      "fit must be an object returned by fpgs().",
      call. = FALSE
    )
  }

  type <- match.arg(type)

  if (is.null(fit$method) ||
      !fit$method %in% c(
        "linear_regression",
        "logistic_regression"
      )) {

    stop(
      "Analytic variance is only available for FPGS fits estimated ",
      "using method = 'linear_regression' or 'logistic_regression'.",
      call. = FALSE
    )
  }

  outcome_type <- fit$outcome_type

  if (is.null(outcome_type)) {
    stop(
      "Outcome type is not stored in the fpgs object.",
      call. = FALSE
    )
  }

  if (!outcome_type %in% c("continuous", "binary")) {
    stop(
      "outcome_type must be 'continuous' or 'binary'.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------
  # Check consistency between method and outcome type
  # ------------------------------------------------------------

  if (fit$method == "linear_regression" &&
      outcome_type != "continuous") {

    stop(
      "method = 'linear_regression' requires a continuous outcome.",
      call. = FALSE
    )
  }

  if (fit$method == "logistic_regression" &&
      outcome_type != "binary") {

    stop(
      "method = 'logistic_regression' requires a binary outcome.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------
  # Extract data
  # ------------------------------------------------------------

  data <- fit$data
  Y <- data[[fit$outcome]]
  Tr <- data[[fit$treatment]]

  # ------------------------------------------------------------
  # First-stage design matrix
  # ------------------------------------------------------------

  X <- stats::model.matrix(
    stats::reformulate(fit$covariates),
    data = data
  )

  K <- ncol(X)
  X_names <- colnames(X)

  # ------------------------------------------------------------
  # Fit the two outcome models used to construct the FPGS
  # ------------------------------------------------------------

  if (outcome_type == "binary") {

    fit1 <- stats::glm.fit(
      x = X[Tr == 1, , drop = FALSE],
      y = Y[Tr == 1],
      family = stats::binomial(),
      intercept = FALSE
    )

    fit0 <- stats::glm.fit(
      x = X[Tr == 0, , drop = FALSE],
      y = Y[Tr == 0],
      family = stats::binomial(),
      intercept = FALSE
    )

  } else {

    fit1 <- stats::lm.fit(
      x = X[Tr == 1, , drop = FALSE],
      y = Y[Tr == 1]
    )

    fit0 <- stats::lm.fit(
      x = X[Tr == 0, , drop = FALSE],
      y = Y[Tr == 0]
    )
  }

  beta1_0 <- fit1$coefficients
  beta0_0 <- fit0$coefficients

  if (anyNA(beta1_0) || anyNA(beta0_0)) {
    stop(
      "The first-stage outcome models contain unidentified coefficients.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------
  # Construct the two FPGS components
  # ------------------------------------------------------------

  eta1_0 <- drop(
    X %*% beta1_0
  )

  eta0_0 <- drop(
    X %*% beta0_0
  )

  if (outcome_type == "binary") {

    m1_0 <- stats::plogis(
      eta1_0
    )

    m0_0 <- stats::plogis(
      eta0_0
    )

  } else {

    m1_0 <- eta1_0
    m0_0 <- eta0_0
  }

  # ------------------------------------------------------------
  # Augmented FPGS design matrix
  #
  # Z = (1, mu0_hat, mu1_hat)
  # ------------------------------------------------------------

  Z <- cbind(
    `(Intercept)` = 1,
    mu0_hat = m0_0,
    mu1_hat = m1_0
  )

  pg <- ncol(Z)

  # ------------------------------------------------------------
  # FPGS-based propensity score model
  # ------------------------------------------------------------

  ps_fit <- stats::glm.fit(
    x = Z,
    y = Tr,
    family = stats::binomial(),
    intercept = FALSE
  )

  gamma0 <- ps_fit$coefficients

  if (anyNA(gamma0)) {
    stop(
      "The FPGS propensity score model contains unidentified coefficients.",
      call. = FALSE
    )
  }

  e_model0 <- stats::plogis(
    drop(Z %*% gamma0)
  )

  # Truncated propensity scores used by the weighting estimators
  e0 <- pmin(
    pmax(e_model0, 1e-6),
    1 - 1e-6
  )

  # ------------------------------------------------------------
  # Starting values for IPW
  # ------------------------------------------------------------

  if (type == "ipw") {

    mu1_0 <- mean(
      Tr * Y / e0
    )

    mu0_0 <- mean(
      (1 - Tr) * Y / (1 - e0)
    )

    tau_0 <- mu1_0 - mu0_0

    theta0 <- c(
      beta1_0,
      beta0_0,
      gamma0,
      mu1_0,
      mu0_0,
      tau_0
    )
  }

  # ------------------------------------------------------------
  # Starting values for NIPW
  # ------------------------------------------------------------

  if (type == "nipw") {

    a1_0 <- mean(
      Tr * Y / e0
    )

    b1_0 <- mean(
      Tr / e0
    )

    a0_0 <- mean(
      (1 - Tr) * Y / (1 - e0)
    )

    b0_0 <- mean(
      (1 - Tr) / (1 - e0)
    )

    mu1_0 <- a1_0 / b1_0
    mu0_0 <- a0_0 / b0_0

    tau_0 <- mu1_0 - mu0_0

    theta0 <- c(
      beta1_0,
      beta0_0,
      gamma0,
      a1_0,
      b1_0,
      a0_0,
      b0_0,
      mu1_0,
      mu0_0,
      tau_0
    )
  }

  # ------------------------------------------------------------
  # Starting values for AIPW
  # ------------------------------------------------------------

  if (type == "aipw") {

    # ----------------------------------------------------------
    # Additional outcome models based on the estimated FPGS
    # ----------------------------------------------------------

    if (outcome_type == "binary") {

      outcome_fit1 <- stats::glm.fit(
        x = Z[Tr == 1, , drop = FALSE],
        y = Y[Tr == 1],
        family = stats::binomial(),
        intercept = FALSE
      )

      outcome_fit0 <- stats::glm.fit(
        x = Z[Tr == 0, , drop = FALSE],
        y = Y[Tr == 0],
        family = stats::binomial(),
        intercept = FALSE
      )

    } else {

      outcome_fit1 <- stats::lm.fit(
        x = Z[Tr == 1, , drop = FALSE],
        y = Y[Tr == 1]
      )

      outcome_fit0 <- stats::lm.fit(
        x = Z[Tr == 0, , drop = FALSE],
        y = Y[Tr == 0]
      )
    }

    alpha1_0 <- outcome_fit1$coefficients
    alpha0_0 <- outcome_fit0$coefficients

    if (anyNA(alpha1_0) || anyNA(alpha0_0)) {
      stop(
        "The FPGS-based outcome models contain unidentified coefficients.",
        call. = FALSE
      )
    }

    outcome_eta1_0 <- drop(
      Z %*% alpha1_0
    )

    outcome_eta0_0 <- drop(
      Z %*% alpha0_0
    )

    if (outcome_type == "binary") {

      outcome_mean1_0 <- stats::plogis(
        outcome_eta1_0
      )

      outcome_mean0_0 <- stats::plogis(
        outcome_eta0_0
      )

    } else {

      outcome_mean1_0 <- outcome_eta1_0
      outcome_mean0_0 <- outcome_eta0_0
    }

    mu1_0 <- mean(
      outcome_mean1_0 +
        Tr * (Y - outcome_mean1_0) / e0
    )

    mu0_0 <- mean(
      outcome_mean0_0 +
        (1 - Tr) * (Y - outcome_mean0_0) / (1 - e0)
    )

    tau_0 <- mu1_0 - mu0_0

    theta0 <- c(
      beta1_0,
      beta0_0,
      alpha1_0,
      alpha0_0,
      gamma0,
      mu1_0,
      mu0_0,
      tau_0
    )
  }

  # ------------------------------------------------------------
  # Data passed to geex
  # ------------------------------------------------------------

  units <- data.frame(
    Y = Y,
    Tr = Tr,
    as.data.frame(
      X,
      check.names = FALSE
    ),
    check.names = FALSE
  )

  # ------------------------------------------------------------
  # Stacked estimating equations
  # ------------------------------------------------------------

  estFUN <- function(data) {

    function(theta) {

      beta1 <- theta[
        seq_len(K)
      ]

      beta0 <- theta[
        K + seq_len(K)
      ]

      x <- as.numeric(
        unlist(
          data[X_names],
          use.names = FALSE
        )
      )

      # --------------------------------------------------------
      # First-stage outcome models
      # --------------------------------------------------------

      eta1 <- sum(
        x * beta1
      )

      eta0 <- sum(
        x * beta0
      )

      if (outcome_type == "binary") {

        m1 <- stats::plogis(
          eta1
        )

        m0 <- stats::plogis(
          eta0
        )

      } else {

        m1 <- eta1
        m0 <- eta0
      }

      z <- c(
        1,
        m0,
        m1
      )

      # --------------------------------------------------------
      # First-stage outcome-model estimating equations
      # --------------------------------------------------------

      psi_beta1 <-
        data$Tr *
        (data$Y - m1) *
        x

      psi_beta0 <-
        (1 - data$Tr) *
        (data$Y - m0) *
        x

      # ========================================================
      # IPW
      # ========================================================

      if (type == "ipw") {

        gamma <- theta[
          2 * K + seq_len(pg)
        ]

        mu1 <- theta[
          2 * K + pg + 1
        ]

        mu0 <- theta[
          2 * K + pg + 2
        ]

        tau <- theta[
          2 * K + pg + 3
        ]

        e_model <- stats::plogis(
          sum(z * gamma)
        )

        e <- min(
          max(e_model, 1e-6),
          1 - 1e-6
        )

        psi_gamma <-
          (data$Tr - e_model) *
          z

        psi_mu1 <-
          mu1 -
          data$Tr * data$Y / e

        psi_mu0 <-
          mu0 -
          (1 - data$Tr) * data$Y / (1 - e)

        psi_tau <-
          tau -
          (mu1 - mu0)

        return(
          c(
            psi_beta1,
            psi_beta0,
            psi_gamma,
            psi_mu1,
            psi_mu0,
            psi_tau
          )
        )
      }

      # ========================================================
      # NIPW
      # ========================================================

      if (type == "nipw") {

        gamma <- theta[
          2 * K + seq_len(pg)
        ]

        a1 <- theta[
          2 * K + pg + 1
        ]

        b1 <- theta[
          2 * K + pg + 2
        ]

        a0 <- theta[
          2 * K + pg + 3
        ]

        b0 <- theta[
          2 * K + pg + 4
        ]

        mu1 <- theta[
          2 * K + pg + 5
        ]

        mu0 <- theta[
          2 * K + pg + 6
        ]

        tau <- theta[
          2 * K + pg + 7
        ]

        e_model <- stats::plogis(
          sum(z * gamma)
        )

        e <- min(
          max(e_model, 1e-6),
          1 - 1e-6
        )

        psi_gamma <-
          (data$Tr - e_model) *
          z

        psi_a1 <-
          a1 -
          data$Tr * data$Y / e

        psi_b1 <-
          b1 -
          data$Tr / e

        psi_a0 <-
          a0 -
          (1 - data$Tr) * data$Y / (1 - e)

        psi_b0 <-
          b0 -
          (1 - data$Tr) / (1 - e)

        psi_mu1 <-
          mu1 -
          a1 / b1

        psi_mu0 <-
          mu0 -
          a0 / b0

        psi_tau <-
          tau -
          (mu1 - mu0)

        return(
          c(
            psi_beta1,
            psi_beta0,
            psi_gamma,
            psi_a1,
            psi_b1,
            psi_a0,
            psi_b0,
            psi_mu1,
            psi_mu0,
            psi_tau
          )
        )
      }

      # ========================================================
      # AIPW
      # ========================================================

      if (type == "aipw") {

        alpha1 <- theta[
          2 * K + seq_len(pg)
        ]

        alpha0 <- theta[
          2 * K + pg + seq_len(pg)
        ]

        gamma <- theta[
          2 * K + 2 * pg + seq_len(pg)
        ]

        mu1 <- theta[
          2 * K + 3 * pg + 1
        ]

        mu0 <- theta[
          2 * K + 3 * pg + 2
        ]

        tau <- theta[
          2 * K + 3 * pg + 3
        ]

        # ------------------------------------------------------
        # Second-stage FPGS-based outcome means
        # ------------------------------------------------------

        outcome_eta1 <- sum(
          z * alpha1
        )

        outcome_eta0 <- sum(
          z * alpha0
        )

        if (outcome_type == "binary") {

          outcome_mean1 <- stats::plogis(
            outcome_eta1
          )

          outcome_mean0 <- stats::plogis(
            outcome_eta0
          )

        } else {

          outcome_mean1 <- outcome_eta1
          outcome_mean0 <- outcome_eta0
        }

        # ------------------------------------------------------
        # FPGS-based propensity score
        # ------------------------------------------------------

        e_model <- stats::plogis(
          sum(z * gamma)
        )

        e <- min(
          max(e_model, 1e-6),
          1 - 1e-6
        )

        # ------------------------------------------------------
        # Second-stage outcome-model estimating equations
        # ------------------------------------------------------

        psi_alpha1 <-
          data$Tr *
          (data$Y - outcome_mean1) *
          z

        psi_alpha0 <-
          (1 - data$Tr) *
          (data$Y - outcome_mean0) *
          z

        # ------------------------------------------------------
        # Propensity-score estimating equation
        # ------------------------------------------------------

        psi_gamma <-
          (data$Tr - e_model) *
          z

        # ------------------------------------------------------
        # AIPW estimating equations
        # ------------------------------------------------------

        psi_mu1 <-
          mu1 -
          (
            outcome_mean1 +
              data$Tr *
              (data$Y - outcome_mean1) / e
          )

        psi_mu0 <-
          mu0 -
          (
            outcome_mean0 +
              (1 - data$Tr) *
              (data$Y - outcome_mean0) / (1 - e)
          )

        psi_tau <-
          tau -
          (mu1 - mu0)

        return(
          c(
            psi_beta1,
            psi_beta0,
            psi_alpha1,
            psi_alpha0,
            psi_gamma,
            psi_mu1,
            psi_mu0,
            psi_tau
          )
        )
      }
    }
  }

  # ------------------------------------------------------------
  # Fit the stacked M-estimation system
  # ------------------------------------------------------------

  geex_fit <- geex::m_estimate(
    estFUN = estFUN,
    data = units,
    roots = theta0,
    compute_roots = FALSE
  )

  # ------------------------------------------------------------
  # Sandwich covariance matrix
  # ------------------------------------------------------------

  V <- geex::vcov(
    geex_fit
  )

  # ------------------------------------------------------------
  # Location of ATE parameter
  # ------------------------------------------------------------

  if (type == "ipw") {

    tau_index <-
      2 * K + pg + 3
  }

  if (type == "nipw") {

    tau_index <-
      2 * K + pg + 7
  }

  if (type == "aipw") {

    tau_index <-
      2 * K + 3 * pg + 3
  }

  variance <- V[
    tau_index,
    tau_index
  ]

  se <- sqrt(
    variance
  )

  # ------------------------------------------------------------
  # Output
  # ------------------------------------------------------------

  list(
    estimate = unname(
      theta0[tau_index]
    ),
    variance = unname(
      variance
    ),
    se = unname(
      se
    )
  )
}
