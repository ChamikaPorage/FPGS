#' Analytic Variance for the FPGS Regression Imputation Estimator
#'
#' Computes the sandwich variance and standard error for the two-stage
#' FPGS regression imputation estimator using stacked estimating equations.
#'
#' The stacked M-estimation system accounts for estimation of:
#' \enumerate{
#'   \item the two treatment-specific regression models used to estimate
#'   the FPGS components;
#'   \item the two treatment-specific second-stage regression models using
#'   the estimated FPGS components; and
#'   \item the regression-imputation average treatment effect.
#' }
#'
#' For continuous outcomes, linear regression models are used in both
#' stages. For binary outcomes, logistic regression models are used in
#' both stages.
#'
#' @param fit An object of class \code{"fpgs"} returned by \code{fpgs()}
#'   with \code{method = "linear_regression"} or
#'   \code{method = "logistic_regression"}.
#'
#' @return A list containing:
#' \describe{
#'   \item{estimate}{The two-stage regression imputation estimate of the ATE.}
#'   \item{variance}{Estimated sandwich variance of the ATE estimator.}
#'   \item{se}{Estimated standard error of the ATE estimator.}
#' }
#'
#' @details
#' Let the first-stage FPGS components be estimated from treatment-specific
#' outcome models. In the second stage, the outcome is regressed on the two
#' estimated FPGS components separately within the treated and untreated
#' groups. The ATE is obtained by averaging the difference between the two
#' second-stage predicted outcomes.
#'
#' The variance estimator is based on a stacked M-estimation system containing
#' estimating equations for both first-stage outcome models, both second-stage
#' outcome models, and the regression-imputation ATE.
#'
#' This analytic variance estimator is available only when the FPGS is
#' estimated using linear regression or logistic regression. Bootstrap
#' variance estimation should be used when the FPGS components are estimated
#' using random forest.
#'
#' @keywords internal
ri_variance <- function(fit) {

  # ------------------------------------------------------------
  # Check input
  # ------------------------------------------------------------

  if (!inherits(fit, "fpgs")) {
    stop(
      "fit must be an object returned by fpgs().",
      call. = FALSE
    )
  }

  if (!fit$method %in% c(
    "linear_regression",
    "logistic_regression"
  )) {
    stop(
      "ri_variance() is available only for FPGS fits estimated ",
      "using method = 'linear_regression' or ",
      "'logistic_regression'.",
      call. = FALSE
    )
  }

  outcome_type <- fit$outcome_type

  if (!outcome_type %in% c("continuous", "binary")) {
    stop(
      "outcome_type must be either 'continuous' or 'binary'.",
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

  dat <- fit$data

  Y <- dat[[fit$outcome]]
  Tr <- dat[[fit$treatment]]

  # ------------------------------------------------------------
  # First-stage design matrix
  # ------------------------------------------------------------

  X <- stats::model.matrix(
    stats::reformulate(fit$covariates),
    data = dat
  )

  K <- ncol(X)
  Xnames <- colnames(X)

  # ------------------------------------------------------------
  # STEP 1:
  # Fit the two regression models used to estimate the FPGS
  # ------------------------------------------------------------

  if (outcome_type == "continuous") {

    first1 <- stats::lm.fit(
      x = X[Tr == 1, , drop = FALSE],
      y = Y[Tr == 1]
    )

    first0 <- stats::lm.fit(
      x = X[Tr == 0, , drop = FALSE],
      y = Y[Tr == 0]
    )

  } else {

    first1 <- stats::glm.fit(
      x = X[Tr == 1, , drop = FALSE],
      y = Y[Tr == 1],
      family = stats::binomial(),
      intercept = FALSE
    )

    first0 <- stats::glm.fit(
      x = X[Tr == 0, , drop = FALSE],
      y = Y[Tr == 0],
      family = stats::binomial(),
      intercept = FALSE
    )
  }

  beta1_0 <- first1$coefficients
  beta0_0 <- first0$coefficients

  if (anyNA(beta1_0) || anyNA(beta0_0)) {
    stop(
      "The first-stage regression models contain unidentified coefficients.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------
  # First-stage FPGS components
  # ------------------------------------------------------------

  eta1_first <- drop(
    X %*% beta1_0
  )

  eta0_first <- drop(
    X %*% beta0_0
  )

  if (outcome_type == "binary") {

    mu1_0 <- stats::plogis(
      eta1_first
    )

    mu0_0 <- stats::plogis(
      eta0_first
    )

  } else {

    mu1_0 <- eta1_first
    mu0_0 <- eta0_first
  }

  # ------------------------------------------------------------
  # STEP 2:
  # Second-stage design matrix
  #
  # Z = (1, mu0_hat, mu1_hat)
  # ------------------------------------------------------------

  Z <- cbind(
    `(Intercept)` = 1,
    mu0_hat = mu0_0,
    mu1_hat = mu1_0
  )

  Q <- ncol(Z)

  # ------------------------------------------------------------
  # Fit second-stage treatment-specific models
  # ------------------------------------------------------------

  if (outcome_type == "continuous") {

    second1 <- stats::lm.fit(
      x = Z[Tr == 1, , drop = FALSE],
      y = Y[Tr == 1]
    )

    second0 <- stats::lm.fit(
      x = Z[Tr == 0, , drop = FALSE],
      y = Y[Tr == 0]
    )

  } else {

    second1 <- stats::glm.fit(
      x = Z[Tr == 1, , drop = FALSE],
      y = Y[Tr == 1],
      family = stats::binomial(),
      intercept = FALSE
    )

    second0 <- stats::glm.fit(
      x = Z[Tr == 0, , drop = FALSE],
      y = Y[Tr == 0],
      family = stats::binomial(),
      intercept = FALSE
    )
  }

  gamma1_0 <- second1$coefficients
  gamma0_0 <- second0$coefficients

  if (anyNA(gamma1_0) || anyNA(gamma0_0)) {
    stop(
      "The second-stage regression models contain unidentified coefficients.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------
  # Second-stage predicted potential outcomes
  # ------------------------------------------------------------

  eta1_second <- drop(
    Z %*% gamma1_0
  )

  eta0_second <- drop(
    Z %*% gamma0_0
  )

  if (outcome_type == "binary") {

    m1_0 <- stats::plogis(
      eta1_second
    )

    m0_0 <- stats::plogis(
      eta0_second
    )

  } else {

    m1_0 <- eta1_second
    m0_0 <- eta0_second
  }

  # ------------------------------------------------------------
  # Initial ATE
  # ------------------------------------------------------------

  tau_0 <- mean(
    m1_0 - m0_0
  )

  # ------------------------------------------------------------
  # Parameter vector
  #
  # theta =
  # (beta1, beta0, gamma1, gamma0, tau)
  # ------------------------------------------------------------

  theta0 <- c(
    beta1_0,
    beta0_0,
    gamma1_0,
    gamma0_0,
    tau_0
  )

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
  # Parameter indices
  # ------------------------------------------------------------

  idx_beta1 <- seq_len(K)

  idx_beta0 <-
    K + seq_len(K)

  idx_gamma1 <-
    2 * K + seq_len(Q)

  idx_gamma0 <-
    2 * K + Q + seq_len(Q)

  idx_tau <-
    2 * K + 2 * Q + 1

  # ------------------------------------------------------------
  # Stacked estimating equations
  # ------------------------------------------------------------

  estFUN <- function(data) {

    function(theta) {

      # --------------------------------------------------------
      # Extract parameters
      # --------------------------------------------------------

      beta1 <- theta[idx_beta1]
      beta0 <- theta[idx_beta0]

      gamma1 <- theta[idx_gamma1]
      gamma0 <- theta[idx_gamma0]

      tau <- theta[idx_tau]

      # --------------------------------------------------------
      # Original covariate vector
      # --------------------------------------------------------

      x <- as.numeric(
        unlist(
          data[Xnames],
          use.names = FALSE
        )
      )

      # --------------------------------------------------------
      # STEP 1:
      # FPGS components as functions of beta
      # --------------------------------------------------------

      eta1_first_i <- sum(
        x * beta1
      )

      eta0_first_i <- sum(
        x * beta0
      )

      if (outcome_type == "binary") {

        mu1_i <- stats::plogis(
          eta1_first_i
        )

        mu0_i <- stats::plogis(
          eta0_first_i
        )

      } else {

        mu1_i <- eta1_first_i
        mu0_i <- eta0_first_i
      }

      # --------------------------------------------------------
      # First-stage estimating equations
      # --------------------------------------------------------

      psi_beta1 <-
        data$Tr *
        (data$Y - mu1_i) *
        x

      psi_beta0 <-
        (1 - data$Tr) *
        (data$Y - mu0_i) *
        x

      # --------------------------------------------------------
      # STEP 2:
      # Second-stage covariate vector
      # --------------------------------------------------------

      z <- c(
        1,
        mu0_i,
        mu1_i
      )

      eta1_second_i <- sum(
        z * gamma1
      )

      eta0_second_i <- sum(
        z * gamma0
      )

      if (outcome_type == "binary") {

        m1_i <- stats::plogis(
          eta1_second_i
        )

        m0_i <- stats::plogis(
          eta0_second_i
        )

      } else {

        m1_i <- eta1_second_i
        m0_i <- eta0_second_i
      }

      # --------------------------------------------------------
      # Second-stage estimating equations
      # --------------------------------------------------------

      psi_gamma1 <-
        data$Tr *
        (data$Y - m1_i) *
        z

      psi_gamma0 <-
        (1 - data$Tr) *
        (data$Y - m0_i) *
        z

      # --------------------------------------------------------
      # ATE estimating equation
      # --------------------------------------------------------

      psi_tau <-
        tau -
        (m1_i - m0_i)

      # --------------------------------------------------------
      # Full stacked estimating-equation vector
      # --------------------------------------------------------

      c(
        psi_beta1,
        psi_beta0,
        psi_gamma1,
        psi_gamma0,
        psi_tau
      )
    }
  }

  # ------------------------------------------------------------
  # M-estimation
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

  variance <- V[
    idx_tau,
    idx_tau
  ]

  se <- sqrt(
    variance
  )

  # ------------------------------------------------------------
  # Output
  # ------------------------------------------------------------

  list(
    estimate = unname(tau_0),
    variance = unname(variance),
    se = unname(se)
  )
}
