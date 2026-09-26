#' Estimate Continuous FPGS Using Linear Regression
#'
#' Estimates the two components of the Full Prognostic Score (FPGS)
#' for a continuous outcome using separate linear regression models in
#' the treated and untreated groups.
#'
#' @param data A data frame containing the observed data.
#' @param outcome Character string giving the name of the continuous outcome
#'   variable.
#' @param treatment Character string giving the name of the treatment
#'   variable.
#' @param covariates Character vector giving the names of the pretreatment
#'   covariates. If \code{NULL}, all variables other than the outcome and
#'   treatment are used.
#'
#' @return An object of class \code{"fpgs"} containing the estimated FPGS
#'   components, \code{mu0_hat} and \code{mu1_hat}, together with
#'   information required for treatment-effect estimation.
#'
#' @keywords internal
fpgs_continuous_linear_regression <- function(data,
                                              outcome,
                                              treatment,
                                              covariates = NULL) {

  dat <- as.data.frame(data)

  if (!outcome %in% names(dat)) {
    stop("Outcome variable not found in data.")
  }

  if (!treatment %in% names(dat)) {
    stop("Treatment variable not found in data.")
  }

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

  if (!all(tr %in% c(0, 1))) {
    stop("Treatment must be coded as 0 and 1.")
  }

  if (!all(c(0, 1) %in% tr)) {
    stop("Both treatment groups must be present.")
  }

  form <- stats::reformulate(
    covariates,
    response = outcome
  )

  fit0 <- stats::lm(
    form,
    data = dat[tr == 0, , drop = FALSE]
  )

  fit1 <- stats::lm(
    form,
    data = dat[tr == 1, , drop = FALSE]
  )

  mu0_hat <- stats::predict(
    fit0,
    newdata = dat
  )

  mu1_hat <- stats::predict(
    fit1,
    newdata = dat
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
    outcome_type = "continuous",
    method = "linear_regression"
  )

  class(out) <- "fpgs"

  out
}
