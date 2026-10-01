#' FPGS Regression Imputation Estimator
#'
#' Estimates the average treatment effect (ATE) using regression imputation
#' based on the estimated Full Prognostic Score (FPGS).
#'
#' For FPGS components estimated using linear or logistic regression,
#' treatment-specific regression models are fitted in the second stage.
#' For continuous outcomes, linear regression is fitted using the estimated
#' FPGS components as predictors. For binary outcomes with logistic-regression
#' FPGS estimation, logistic regression is fitted using the corresponding
#' linear predictors (logit-transformed FPGS components).
#'
#' For FPGS components estimated using random forest, random forest with
#' cross-fitting is used in the second stage.
#'
#' For externally estimated FPGS components, treatment-specific regression
#' models are fitted using linear regression for continuous outcomes and
#' logistic regression for binary outcomes.
#'
#' @param fit An object of class \code{"fpgs"} returned by \code{fpgs()},
#'   or an object of class \code{"fpgs_external"} containing externally
#'   estimated FPGS components.
#'
#' @return An object of class \code{"fpgs_ri"} containing the ATE estimate
#'   and predicted potential outcomes.
#'
#' @export
fpgs_ri <- function(fit) {

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
  # Check outcome type
  # ----------------------------------------------------------

  if (is.null(fit$outcome_type) ||
      !fit$outcome_type %in% c("continuous", "binary")) {

    stop(
      "outcome_type must be either 'continuous' or 'binary'.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Extract data
  # ----------------------------------------------------------

  dat <- fit$data
  Y <- dat[[fit$outcome]]
  Tr <- dat[[fit$treatment]]

  fpgs_dat <- data.frame(
    Y = Y,
    Tr = Tr,
    mu0_hat = fit$mu0_hat,
    mu1_hat = fit$mu1_hat
  )

  # ----------------------------------------------------------
  # Add linear predictors for logistic-regression FPGS
  # ----------------------------------------------------------

  if (fit$method == "logistic_regression") {

    if (is.null(fit$eta0_hat) || is.null(fit$eta1_hat)) {
      stop(
        "Linear predictors eta0_hat and eta1_hat are required ",
        "for logistic-regression FPGS.",
        call. = FALSE
      )
    }

    fpgs_dat$eta0_hat <- fit$eta0_hat
    fpgs_dat$eta1_hat <- fit$eta1_hat
  }

  # ----------------------------------------------------------
  # Check FPGS values
  # ----------------------------------------------------------

  if (anyNA(fpgs_dat)) {
    stop(
      "Outcome, treatment, and estimated FPGS components must not ",
      "contain missing values.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Regression-based or externally supplied FPGS
  # ----------------------------------------------------------

  if (fit$method %in% c(
    "linear_regression",
    "logistic_regression",
    "external"
  )) {

    # --------------------------------------------------------
    # Continuous outcome
    # --------------------------------------------------------

    if (fit$outcome_type == "continuous") {

      model0 <- stats::lm(
        Y ~ mu0_hat + mu1_hat,
        data = fpgs_dat[Tr == 0, , drop = FALSE]
      )

      model1 <- stats::lm(
        Y ~ mu0_hat + mu1_hat,
        data = fpgs_dat[Tr == 1, , drop = FALSE]
      )

      pred0 <- stats::predict(
        model0,
        newdata = fpgs_dat
      )

      pred1 <- stats::predict(
        model1,
        newdata = fpgs_dat
      )

      second_stage <- "linear regression"
    }

    # --------------------------------------------------------
    # Binary outcome
    # --------------------------------------------------------

    if (fit$outcome_type == "binary") {

      # ------------------------------------------------------
      # Logistic-regression FPGS:
      # use linear predictors (logit-transformed FPGS)
      # ------------------------------------------------------

      if (fit$method == "logistic_regression") {

        model0 <- stats::glm(
          Y ~ eta0_hat + eta1_hat,
          data = fpgs_dat[Tr == 0, , drop = FALSE],
          family = stats::binomial()
        )

        model1 <- stats::glm(
          Y ~ eta0_hat + eta1_hat,
          data = fpgs_dat[Tr == 1, , drop = FALSE],
          family = stats::binomial()
        )
      }

      # ------------------------------------------------------
      # Externally supplied FPGS:
      # use supplied FPGS components
      # ------------------------------------------------------

      if (fit$method == "external") {

        model0 <- stats::glm(
          Y ~ mu0_hat + mu1_hat,
          data = fpgs_dat[Tr == 0, , drop = FALSE],
          family = stats::binomial()
        )

        model1 <- stats::glm(
          Y ~ mu0_hat + mu1_hat,
          data = fpgs_dat[Tr == 1, , drop = FALSE],
          family = stats::binomial()
        )
      }

      pred0 <- stats::predict(
        model0,
        newdata = fpgs_dat,
        type = "response"
      )

      pred1 <- stats::predict(
        model1,
        newdata = fpgs_dat,
        type = "response"
      )

      second_stage <- "logistic regression"
    }

    # --------------------------------------------------------
    # Average treatment effect
    # --------------------------------------------------------

    ate <- mean(pred1 - pred0)

    out <- list(
      estimate = ate,
      pred0 = pred0,
      pred1 = pred1,
      model0 = model0,
      model1 = model1,
      second_stage = second_stage,
      method = "FPGS regression imputation"
    )
  }

  # ----------------------------------------------------------
  # Random-forest FPGS
  # ----------------------------------------------------------

  if (fit$method == "random_forest") {

    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop(
        "Package 'ranger' is required.",
        call. = FALSE
      )
    }

    if (!requireNamespace("caret", quietly = TRUE)) {
      stop(
        "Package 'caret' is required.",
        call. = FALSE
      )
    }

    folds <- fit$folds

    if (is.null(folds)) {
      folds <- 5
    }

    num.trees <- fit$num.trees

    if (is.null(num.trees)) {
      num.trees <- 500
    }

    n <- nrow(fpgs_dat)

    fold_id <- caret::createFolds(
      Y,
      k = folds,
      list = TRUE
    )

    pred0 <- rep(NA_real_, n)
    pred1 <- rep(NA_real_, n)

    for (k in seq_len(folds)) {

      train_index <- unlist(fold_id[-k])
      test_index <- fold_id[[k]]

      train_data <- fpgs_dat[
        train_index,
        ,
        drop = FALSE
      ]

      test_data <- fpgs_dat[
        test_index,
        ,
        drop = FALSE
      ]

      train0 <- train_data[
        train_data$Tr == 0,
        ,
        drop = FALSE
      ]

      train1 <- train_data[
        train_data$Tr == 1,
        ,
        drop = FALSE
      ]

      # ------------------------------------------------------
      # Continuous outcome
      # ------------------------------------------------------

      if (fit$outcome_type == "continuous") {

        model0 <- ranger::ranger(
          Y ~ mu0_hat + mu1_hat,
          data = train0,
          num.trees = num.trees,
          mtry = 2
        )

        model1 <- ranger::ranger(
          Y ~ mu0_hat + mu1_hat,
          data = train1,
          num.trees = num.trees,
          mtry = 2
        )

        pred0[test_index] <- stats::predict(
          model0,
          data = test_data
        )$predictions

        pred1[test_index] <- stats::predict(
          model1,
          data = test_data
        )$predictions
      }

      # ------------------------------------------------------
      # Binary outcome
      # ------------------------------------------------------

      if (fit$outcome_type == "binary") {

        train0$Y <- factor(
          train0$Y,
          levels = c(0, 1)
        )

        train1$Y <- factor(
          train1$Y,
          levels = c(0, 1)
        )

        model0 <- ranger::ranger(
          Y ~ mu0_hat + mu1_hat,
          data = train0,
          probability = TRUE,
          num.trees = num.trees,
          mtry = 2
        )

        model1 <- ranger::ranger(
          Y ~ mu0_hat + mu1_hat,
          data = train1,
          probability = TRUE,
          num.trees = num.trees,
          mtry = 2
        )

        prediction0 <- stats::predict(
          model0,
          data = test_data
        )$predictions

        prediction1 <- stats::predict(
          model1,
          data = test_data
        )$predictions

        pred0[test_index] <- prediction0[, "1"]
        pred1[test_index] <- prediction1[, "1"]
      }
    }

    # --------------------------------------------------------
    # Average treatment effect
    # --------------------------------------------------------

    ate <- mean(pred1 - pred0)

    out <- list(
      estimate = ate,
      pred0 = pred0,
      pred1 = pred1,
      folds = folds,
      num.trees = num.trees,
      second_stage = "random forest",
      method = "FPGS regression imputation"
    )
  }

  # ----------------------------------------------------------
  # Check FPGS estimation method
  # ----------------------------------------------------------

  if (!fit$method %in% c(
    "linear_regression",
    "logistic_regression",
    "random_forest",
    "external"
  )) {

    stop(
      "Unknown FPGS estimation method.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Output
  # ----------------------------------------------------------

  class(out) <- "fpgs_ri"

  out
}
