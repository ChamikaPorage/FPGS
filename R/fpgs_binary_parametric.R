
#' Estimate Binary FPGS Using Logistic Regression
#'
#' Estimates the two components of the Full Prognostic Score (FPGS)
#' for a binary outcome using separate logistic regression models in
#' the treated and untreated groups.
#'
#' @param data A data frame containing the observed data.
#' @param outcome Character string giving the name of the binary outcome
#'   variable.
#' @param treatment Character string giving the name of the treatment
#'   variable.
#' @param covariates Character vector giving the names of the pretreatment
#'   covariates. If \code{NULL}, all variables other than the outcome and
#'   treatment are used.
#'
#' @return An object of class \code{"fpgs"} containing the estimated FPGS
#'   components and information required for treatment-effect estimation.
#'
#' @keywords internal
fpgs_binary_logistic_regression <- function(data,
                                            outcome,
                                            treatment,
                                            covariates = NULL) {

  dat <- as.data.frame(data)

  if (is.null(covariates)) {
    covariates <- setdiff(
      names(dat),
      c(outcome, treatment)
    )
  }

  if (length(covariates) == 0) {
    stop("At least one pretreatment covariate must be supplied.")
  }

  tr <- dat[[treatment]]

  form <- stats::reformulate(
    covariates,
    response = outcome
  )

  fit0 <- stats::glm(
    form,
    data = dat[tr == 0, , drop = FALSE],
    family = stats::binomial()
  )

  fit1 <- stats::glm(
    form,
    data = dat[tr == 1, , drop = FALSE],
    family = stats::binomial()
  )

  mu0_hat <- stats::predict(
    fit0,
    newdata = dat,
    type = "response"
  )

  mu1_hat <- stats::predict(
    fit1,
    newdata = dat,
    type = "response"
  )

  out <- list(
    data = dat,
    outcome = outcome,
    treatment = treatment,
    covariates = covariates,
    mu0_hat = mu0_hat,
    mu1_hat = mu1_hat,
    fpgs = data.frame(
      mu0_hat = mu0_hat,
      mu1_hat = mu1_hat
    ),
    outcome_type = "binary",
    method = "logistic_regression"
  )

  class(out) <- "fpgs"

  out
}
