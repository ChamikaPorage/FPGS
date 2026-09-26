#' Print an FPGS Fit
#'
#' Prints a brief summary of an object returned by \code{fpgs()},
#' including the outcome, treatment, outcome type, FPGS estimation method,
#' and number of covariates.
#'
#' @param x An object of class \code{"fpgs"} returned by \code{fpgs()}.
#' @param ... Additional arguments passed to the print method.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.fpgs <- function(x, ...) {

  cat("\nFull Prognostic Score Model\n")
  cat("----------------------------------\n")

  cat("Outcome:", x$outcome, "\n")
  cat("Treatment:", x$treatment, "\n")

  if (!is.null(x$outcome_type)) {
    cat("Outcome type:", x$outcome_type, "\n")
  }

  if (!is.null(x$method)) {

    method_label <- switch(
      x$method,
      linear_regression = "linear regression",
      logistic_regression = "logistic regression",
      random_forest = "random forest",
      x$method
    )

    cat("FPGS estimation method:", method_label, "\n")
  }

  if (!is.null(x$covariates)) {
    cat("Number of covariates:", length(x$covariates), "\n")
  }

  invisible(x)
}
