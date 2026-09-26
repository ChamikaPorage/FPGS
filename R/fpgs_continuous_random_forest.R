#' Estimate Continuous FPGS Using Random Forest with Cross-Fitting
#'
#' Estimates the two components of the Full Prognostic Score (FPGS)
#' for a continuous outcome using separate random forest models in the
#' treated and untreated groups with cross-fitting.
#'
#' @param data A data frame containing the observed data.
#' @param outcome Character string giving the name of the continuous outcome
#'   variable.
#' @param treatment Character string giving the name of the treatment
#'   variable.
#' @param covariates Character vector giving the names of the pretreatment
#'   covariates. If \code{NULL}, all variables other than the outcome and
#'   treatment are used.
#' @param folds Number of folds used for cross-fitting. Default is 5.
#' @param num.trees Number of trees used in each random forest. Default is 500.
#' @param mtry Number of variables considered at each split. If \code{NULL},
#'   it is set to the square root of the number of covariates.
#' @param min.node.size Minimum terminal node size. Default is 5.
#'
#' @return An object of class \code{"fpgs"} containing the cross-fitted
#'   estimates of the two FPGS components, \code{mu0_hat} and
#'   \code{mu1_hat}, together with information required for
#'   treatment-effect estimation.
#'
#' @keywords internal
fpgs_continuous_random_forest <- function(data,
                                          outcome,
                                          treatment,
                                          covariates = NULL,
                                          folds = 5,
                                          num.trees = 500,
                                          mtry = NULL,
                                          min.node.size = 5) {

  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required.")
  }

  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is required.")
  }

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

  if (is.null(mtry)) {
    mtry <- max(
      1,
      floor(sqrt(length(covariates)))
    )
  }

  Y <- dat[[outcome]]
  Tr <- dat[[treatment]]
  n <- nrow(dat)

  if (!all(Tr %in% c(0, 1))) {
    stop("Treatment must be coded as 0 and 1.")
  }

  if (!all(c(0, 1) %in% Tr)) {
    stop("Both treatment groups must be present.")
  }

  fold_id <- caret::createFolds(
    Y,
    k = folds,
    list = TRUE
  )

  mu0_hat <- rep(NA_real_, n)
  mu1_hat <- rep(NA_real_, n)

  model_data <- dat[
    ,
    c(outcome, treatment, covariates),
    drop = FALSE
  ]

  form <- stats::reformulate(
    covariates,
    response = outcome
  )

  for (k in seq_len(folds)) {

    train_index <- unlist(fold_id[-k])
    valid_index <- fold_id[[k]]

    train_data <- model_data[
      train_index,
      ,
      drop = FALSE
    ]

    valid_data <- model_data[
      valid_index,
      ,
      drop = FALSE
    ]

    train_data_0 <- train_data[
      train_data[[treatment]] == 0,
      ,
      drop = FALSE
    ]

    train_data_1 <- train_data[
      train_data[[treatment]] == 1,
      ,
      drop = FALSE
    ]

    rf0 <- ranger::ranger(
      formula = form,
      data = train_data_0,
      num.trees = num.trees,
      mtry = mtry,
      min.node.size = min.node.size
    )

    rf1 <- ranger::ranger(
      formula = form,
      data = train_data_1,
      num.trees = num.trees,
      mtry = mtry,
      min.node.size = min.node.size
    )

    mu0_hat[valid_index] <- stats::predict(
      rf0,
      data = valid_data
    )$predictions

    mu1_hat[valid_index] <- stats::predict(
      rf1,
      data = valid_data
    )$predictions
  }

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
    method = "random_forest",
    learner = "random forest",
    crossfit = TRUE,
    folds = folds,
    num.trees = num.trees,
    mtry = mtry,
    min.node.size = min.node.size
  )

  class(out) <- "fpgs"

  out
}
