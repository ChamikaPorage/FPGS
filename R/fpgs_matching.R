#' FPGS Matching Estimator
#'
#' Estimates the average treatment effect (ATE) by matching on the
#' two-dimensional estimated Full Prognostic Score (FPGS).
#'
#' @param fit An object of class \code{"fpgs"} returned by \code{fpgs()},
#'   or an object of class \code{"fpgs_external"} containing externally
#'   estimated FPGS components.
#' @param M Number of matches used for each unit. Must be a positive integer.
#'   Default is 1.
#'
#' @return An object of class \code{"fpgs_match"} containing the ATE estimate,
#'   estimated variance and standard error when available, number of matches,
#'   and the matching object returned by \code{Matching::Match()}.
#'
#' @details
#' Matching is performed directly on the two-dimensional estimated FPGS using
#' \code{Matching::Match()}. For FPGS objects returned by \code{fpgs()}, the
#' standard error is obtained from \code{Matching::Match()}, and the variance
#' is calculated as the squared standard error.
#'
#' The FPGS components may also be supplied externally. For externally
#' estimated FPGS components, the matching estimator is available for point
#' estimation, but the variance estimators currently implemented in
#' \code{FPGS} are not available.
#'
#' Matching is not included in the ordinary nonparametric bootstrap
#' implemented by \code{fpgs_bootstrap()}.
#'
#' @export
fpgs_match <- function(fit, M = 1) {

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
  # Check Matching package
  # ----------------------------------------------------------

  if (!requireNamespace("Matching", quietly = TRUE)) {
    stop(
      "Package 'Matching' is required. ",
      "Install it with install.packages('Matching').",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Check number of matches
  # ----------------------------------------------------------

  if (!is.numeric(M) ||
      length(M) != 1 ||
      !is.finite(M) ||
      M < 1 ||
      M != as.integer(M)) {

    stop(
      "M must be a positive integer.",
      call. = FALSE
    )
  }

  M <- as.integer(M)

  # ----------------------------------------------------------
  # Extract data
  # ----------------------------------------------------------

  dat <- fit$data

  Y <- dat[[fit$outcome]]
  Tr <- dat[[fit$treatment]]

  # ----------------------------------------------------------
  # Construct two-dimensional FPGS
  # ----------------------------------------------------------

  X <- cbind(
    mu0_hat = fit$mu0_hat,
    mu1_hat = fit$mu1_hat
  )

  storage.mode(X) <- "double"

  if (anyNA(X)) {
    stop(
      "The estimated FPGS components must not contain missing values.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Matching
  # ----------------------------------------------------------

  match_obj <- Matching::Match(
    Y = Y,
    Tr = Tr,
    X = X,
    estimand = "ATE",
    M = M
  )

  estimate <- as.numeric(match_obj$est)

  # ----------------------------------------------------------
  # Standard error and variance
  # ----------------------------------------------------------

  if (inherits(fit, "fpgs")) {

    se <- unname(match_obj$se.standard)
    variance <- se^2

  } else {

    se <- NA_real_
    variance <- NA_real_
  }

  # ----------------------------------------------------------
  # Output
  # ----------------------------------------------------------

  out <- list(
    estimate = estimate,
    variance = variance,
    se = se,
    M = M,
    match_object = match_obj,
    method = "FPGS matching",
    estimand = "ATE"
  )

  class(out) <- "fpgs_match"

  out
}
