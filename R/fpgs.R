#' Estimate the Full Prognostic Score
#'
#' Estimates the two-dimensional Full Prognostic Score (FPGS) using
#' linear regression, logistic regression, or random forest outcome models
#' for continuous or binary outcomes.
#'
#' @param data A data frame containing the observed data.
#' @param outcome Character string giving the name of the outcome variable.
#' @param treatment Character string giving the name of the treatment variable.
#' @param covariates Character vector giving the names of the pretreatment
#'   covariates. If \code{NULL}, all variables other than the outcome and
#'   treatment are used.
#' @param outcome_type Type of outcome: \code{"continuous"} or \code{"binary"}.
#' @param method FPGS estimation method. Available methods are
#'   \code{"linear_regression"} for continuous outcomes,
#'   \code{"logistic_regression"} for binary outcomes, and
#'   \code{"random_forest"} for either continuous or binary outcomes.
#' @param folds Number of folds used for cross-fitting when
#'   \code{method = "random_forest"}. Default is 5.
#' @param ... Additional arguments passed to the random forest fitting
#'   functions when \code{method = "random_forest"}, such as
#'   \code{num.trees}, \code{mtry}, and \code{min.node.size}.
#'
#' @return An object of class \code{"fpgs"} containing the estimated
#'   FPGS and information required for treatment-effect estimation
#'   and the implemented variance estimators.
#'
#' @export
fpgs <- function(data,
                 outcome,
                 treatment,
                 covariates = NULL,
                 outcome_type = c("continuous", "binary"),
                 method = c(
                   "linear_regression",
                   "logistic_regression",
                   "random_forest"
                 ),
                 folds = 5,
                 ...) {

  outcome_type <- match.arg(outcome_type)
  method <- match.arg(method)

  # Check that the selected method is appropriate for the outcome type
  if (outcome_type == "continuous" &&
      method == "logistic_regression") {
    stop(
      "'logistic_regression' is only available for binary outcomes."
    )
  }

  if (outcome_type == "binary" &&
      method == "linear_regression") {
    stop(
      "'linear_regression' is only available for continuous outcomes."
    )
  }

  # Continuous outcome: linear regression
  if (outcome_type == "continuous" &&
      method == "linear_regression") {

    fit <- fpgs_continuous_linear_regression(
      data,
      outcome,
      treatment,
      covariates
    )
  }

  # Continuous outcome: random forest
  if (outcome_type == "continuous" &&
      method == "random_forest") {

    fit <- fpgs_continuous_random_forest(
      data,
      outcome,
      treatment,
      covariates,
      folds = folds,
      ...
    )
  }

  # Binary outcome: logistic regression
  if (outcome_type == "binary" &&
      method == "logistic_regression") {

    fit <- fpgs_binary_logistic_regression(
      data,
      outcome,
      treatment,
      covariates
    )
  }

  # Binary outcome: random forest
  if (outcome_type == "binary" &&
      method == "random_forest") {

    fit <- fpgs_binary_random_forest(
      data,
      outcome,
      treatment,
      covariates,
      folds = folds,
      ...
    )
  }

  # Store information needed for treatment-effect estimation
  fit$outcome_type <- outcome_type
  fit$method <- method

  fit
}
