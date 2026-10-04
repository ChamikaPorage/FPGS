# ============================================================
# Externally estimated FPGS using XGBoost (simulation example)
# ============================================================

# Load Packages
library(xgboost)
library(caret)
library(FPGS)

# ------------------------------------------------------------
# 1. Generate one simulated dataset
# ------------------------------------------------------------

set.seed(2026)

N <- 1000

# Baseline covariates
X1 <- rnorm(N)
X2 <- rnorm(N)
X3 <- rnorm(N)
X4 <- rnorm(N)
X5 <- rbinom(N, 1, 0.5)
X6 <- rnorm(N)

# Binary effect modifiers
E1 <- rbinom(N, 1, 0.5)
E2 <- rbinom(N, 1, 0.5)
E3 <- rbinom(N, 1, 0.5)

# Treatment assignment
lp_t <- 0.1 * X1 -
  0.1 * X2 +
  1.1 * X3 -
  1.1 * X4 +
  0.4 * X5

p_t <- plogis(lp_t)

Tr <- rbinom(N, 1, p_t)

# Continuous outcome
mu <- -3.85 +
  0.5 * X1 -
  2 * X2 -
  0.5 * X3 +
  2 * X4 +
  X6 -
  E1 -
  2 * E3 +
  5 * Tr +
  Tr * E1 +
  4 * Tr * E2 -
  4 * Tr * E3

Y <- rnorm(N, mean = mu, sd = 1)

# Final dataset
dat <- data.frame(
  Y = Y,
  Tr = Tr,
  X1 = X1,
  X2 = X2,
  X3 = X3,
  X4 = X4,
  X5 = X5,
  X6 = X6,
  E1 = E1,
  E2 = E2,
  E3 = E3
)

# Inspect the data
dim(dat)
head(dat)
table(dat$Tr)

# ------------------------------------------------------------
# 2. True treatment effect
# ------------------------------------------------------------

# Individual treatment effects
true_effect <- 5 + E1 + 4 * E2 - 4 * E3

# Population ATE
true_ate <- 5.5

# Realized sample-average treatment effect
sample_ate <- mean(true_effect)

true_ate
sample_ate

# ------------------------------------------------------------
# 3. Set up 5-fold cross-fitting
# ------------------------------------------------------------

K <- 5

folds <- createFolds(
  dat$Y,
  k = K,
  list = TRUE,
  returnTrain = FALSE
)

# Covariates
x_var <- c(
  "X1", "X2", "X3", "X4", "X5",
  "X6", "E1", "E2", "E3"
)

# Storage for out-of-fold predictions
mu0_xgb <- rep(NA_real_, N)
mu1_xgb <- rep(NA_real_, N)

# ------------------------------------------------------------
# 4. Estimate external FPGS components using XGBoost
# ------------------------------------------------------------

for (k in seq_len(K)) {

  # Held-out fold
  valid_id <- folds[[k]]

  # Training observations
  train_id <- setdiff(seq_len(N), valid_id)

  train_data <- dat[train_id, ]
  valid_data <- dat[valid_id, ]

  # Separate training observations by treatment
  data0_train <- train_data[train_data$Tr == 0, ]
  data1_train <- train_data[train_data$Tr == 1, ]

  # Numeric matrices
  X0_train <- data.matrix(
    data0_train[, x_var, drop = FALSE]
  )

  X1_train <- data.matrix(
    data1_train[, x_var, drop = FALSE]
  )

  X_valid <- data.matrix(
    valid_data[, x_var, drop = FALSE]
  )

  # ----------------------------------------------------------
  # Outcome model under control: T = 0
  # ----------------------------------------------------------

  model0 <- xgboost(
    x = X0_train,
    y = data0_train$Y,
    objective = "reg:squarederror",
    max_depth = 6,
    learning_rate = 0.05,
    nrounds = 300,
    min_child_weight = 5,
    subsample = 0.8,
    colsample_bytree = 0.8
  )

  # ----------------------------------------------------------
  # Outcome model under treatment: T = 1
  # ----------------------------------------------------------

  model1 <- xgboost(
    x = X1_train,
    y = data1_train$Y,
    objective = "reg:squarederror",
    max_depth = 6,
    learning_rate = 0.05,
    nrounds = 500,
    min_child_weight = 5,
    subsample = 0.8,
    colsample_bytree = 0.8
  )

  # ----------------------------------------------------------
  # Out-of-fold predictions
  # ----------------------------------------------------------

  mu0_xgb[valid_id] <- predict(
    model0,
    X_valid
  )

  mu1_xgb[valid_id] <- predict(
    model1,
    X_valid
  )
}

# ------------------------------------------------------------
# 5. Check cross-fitted predictions
# ------------------------------------------------------------

length(mu0_xgb)
length(mu1_xgb)

anyNA(mu0_xgb)
anyNA(mu1_xgb)

summary(mu0_xgb)
summary(mu1_xgb)

# ------------------------------------------------------------
# 6. Estimate ATE using external FPGS components
# ------------------------------------------------------------

ate_xgb <- fpgs_ate(
  mu0 = mu0_xgb,
  mu1 = mu1_xgb,
  data = dat,
  outcome = "Y",
  treatment = "Tr",
  outcome_type = "continuous",
  variance = "none"
)

ate_xgb

# ------------------------------------------------------------
# 7. Compare estimates with the true ATE
# ------------------------------------------------------------

ate_xgb

cat("\nPopulation ATE:", true_ate, "\n")
cat("Sample ATE:", sample_ate, "\n")
