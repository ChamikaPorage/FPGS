test_that("fpgs_ate returns all six estimators", {

  set.seed(2025)

  n <- 500

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

  result <- fpgs_ate(
    fit = fit,
    variance = "none",
    n_bins = 2
  )

  expect_s3_class(result, "data.frame")

  expect_equal(
    names(result),
    c("estimator", "estimate")
  )

  expect_equal(
    result$estimator,
    c(
      "RI",
      "Matching",
      "IPW",
      "NIPW",
      "AIPW",
      "Stratification"
    )
  )

  expect_equal(nrow(result), 6)

  expect_true(all(is.finite(result$estimate)))
})
