test_that("fpgs estimates FPGS for a continuous outcome", {

  set.seed(2025)

  n <- 200

  dat <- data.frame(
    T  = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  dat$Y <- 2 +
    1.5 * dat$T +
    0.5 * dat$X1 -
    0.3 * dat$X2 +
    rnorm(n)

  fit <- fpgs(
    data = dat,
    outcome = "Y",
    treatment = "T",
    covariates = c("X1", "X2"),
    outcome_type = "continuous",
    method = "linear_regression"
  )

  # Correct object class
  expect_s3_class(fit, "fpgs")

  # Correct estimation information
  expect_equal(fit$outcome_type, "continuous")
  expect_equal(fit$method, "linear_regression")

  # One pair of FPGS components for every observation
  expect_length(fit$mu0_hat, n)
  expect_length(fit$mu1_hat, n)

  # Predictions should be available for every observation
  expect_false(anyNA(fit$mu0_hat))
  expect_false(anyNA(fit$mu1_hat))
})

test_that("fpgs estimates FPGS for a binary outcome", {

  set.seed(2026)

  n <- 300

  dat <- data.frame(
    T  = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  lp <- -0.5 +
    0.5 * dat$T +
    0.4 * dat$X1 -
    0.3 * dat$X2

  dat$Y <- rbinom(
    n,
    size = 1,
    prob = plogis(lp)
  )

  fit <- fpgs(
    data = dat,
    outcome = "Y",
    treatment = "T",
    covariates = c("X1", "X2"),
    outcome_type = "binary",
    method = "logistic_regression"
  )

  # Correct object class
  expect_s3_class(fit, "fpgs")

  # Correct estimation information
  expect_equal(fit$outcome_type, "binary")
  expect_equal(fit$method, "logistic_regression")

  # One pair of FPGS components for every observation
  expect_length(fit$mu0_hat, n)
  expect_length(fit$mu1_hat, n)

  # Predictions should be available for every observation
  expect_false(anyNA(fit$mu0_hat))
  expect_false(anyNA(fit$mu1_hat))

  # For a binary outcome, estimated conditional means are probabilities
  expect_true(all(fit$mu0_hat >= 0 & fit$mu0_hat <= 1))
  expect_true(all(fit$mu1_hat >= 0 & fit$mu1_hat <= 1))
})

test_that("fpgs estimates FPGS using random forest with cross-fitting", {

  set.seed(2027)

  n <- 200

  dat <- data.frame(
    T  = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  dat$Y <- 2 +
    1.5 * dat$T +
    0.5 * dat$X1 -
    0.3 * dat$X2 +
    rnorm(n)

  fit <- fpgs(
    data = dat,
    outcome = "Y",
    treatment = "T",
    covariates = c("X1", "X2"),
    outcome_type = "continuous",
    method = "random_forest",
    folds = 3,
    num.trees = 50
  )

  # Correct object class
  expect_s3_class(fit, "fpgs")

  # Correct estimation information
  expect_equal(fit$outcome_type, "continuous")
  expect_equal(fit$method, "random_forest")

  # Cross-fitted predictions for every observation
  expect_length(fit$mu0_hat, n)
  expect_length(fit$mu1_hat, n)

  expect_false(anyNA(fit$mu0_hat))
  expect_false(anyNA(fit$mu1_hat))

  expect_true(all(is.finite(fit$mu0_hat)))
  expect_true(all(is.finite(fit$mu1_hat)))
})

test_that("fpgs rejects logistic regression for continuous outcomes", {

  expect_error(
    fpgs(
      data = dat,
      outcome = "Y",
      treatment = "T",
      covariates = c("X1", "X2"),
      outcome_type = "continuous",
      method = "logistic_regression"
    ),
    "'logistic_regression' is only available for binary outcomes.",
    fixed = TRUE
  )
})


test_that("fpgs rejects linear regression for binary outcomes", {

  set.seed(456)

  n <- 200

  dat_bin <- data.frame(
    T  = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  dat_bin$Y <- rbinom(
    n,
    size = 1,
    prob = plogis(
      -0.5 +
        0.5 * dat_bin$T +
        0.4 * dat_bin$X1 -
        0.3 * dat_bin$X2
    )
  )

  expect_error(
    fpgs(
      data = dat_bin,
      outcome = "Y",
      treatment = "T",
      covariates = c("X1", "X2"),
      outcome_type = "binary",
      method = "linear_regression"
    ),
    "'linear_regression' is only available for continuous outcomes.",
    fixed = TRUE
  )
})
