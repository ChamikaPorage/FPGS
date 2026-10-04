# ============================================================
# FPGS package illustration I
# Continuous outcome: Blood lead
# Treatment: Smoking
# ============================================================

###  Load packages

library(dplyr)
library(FPGS)

### Read data

datat <- read.csv(
  "C:/Users/chapo752/Dropbox/PhD work- Chamika Porage/Second paper/R codes/Emp_ana_new/df_all.csv"
)

###  Recode gender

datat <- datat %>%
  mutate(
    gender = case_when(
      gender == 1 ~ 1,
      gender == 2 ~ 0,
      TRUE ~ NA_real_
    )
  )



###  Covariate matrix


X.matrix <- model.matrix(
  smoking ~ married.living.with.partner +
    birth.country +
    edu +
    race +
    income +
    army.service +
    c.age +
    c.age2 +
    c.family.size +
    gender,
  data = datat
)

# Convert model matrix to data frame
dat.X <- as.data.frame(X.matrix)

# Remove intercept
dat.X <- dat.X[, -1, drop = FALSE]

# Remove incomeOther, as in the previous analysis
dat.X <- dat.X[
  ,
  !names(dat.X) %in% "incomeOther",
  drop = FALSE
]

# Dataset


df_cont <- data.frame(
  smoking = datat$smoking,
  lead = datat$lead,
  dat.X
)

# Inspecting data

dim(df_cont)

names(df_cont)

head(df_cont)

table(df_cont$smoking)

summary(df_cont$lead)

colSums(is.na(df_cont))

# Remove observations with missing values (if present)

df_cont <- na.omit(df_cont)

dim(df_cont)

table(df_cont$smoking)



# Covariates

covariates <- setdiff(
  names(df_cont),
  c("smoking", "lead")
)

length(covariates)

covariates


# ============================================================
# LINEAR REGRESSION FPGS
# ============================================================

#Estimate FPGS using linear regression


fit_cont <- fpgs(
  data = df_cont,
  outcome = "lead",
  treatment = "smoking",
  covariates = covariates,
  outcome_type = "continuous",
  method = "linear_regression"
)

fit_cont

# Check estimated FPGS components


summary(fit_cont$mu0_hat)

summary(fit_cont$mu1_hat)

length(fit_cont$mu0_hat)

length(fit_cont$mu1_hat)

anyNA(fit_cont$mu0_hat)

anyNA(fit_cont$mu1_hat)

### POINT ESTIMATION ###

### Estimate ATE using all six estimators


ate_cont <- fpgs_ate(
  fit = fit_cont,
  variance = "none"
)

ate_cont


#Analytic variance estimation


ate_cont_analytic <- fpgs_ate(
  fit = fit_cont,
  variance = "analytic"
)

ate_cont_analytic



# Checking for analytic results


names(ate_cont_analytic)

all(
  is.finite(ate_cont_analytic$estimate)
)

all(
  ate_cont_analytic$variance >= 0,
  na.rm = TRUE
)

all(
  ate_cont_analytic$se >= 0,
  na.rm = TRUE
)

# Bootstrap variance estimation

set.seed(2026)

ate_cont_boot <- fpgs_ate(
  fit = fit_cont,
  variance = "bootstrap",
  B = 1000,
  conf.level = 0.95
)

ate_cont_boot


# ============================================================
# RANDOM-FOREST FPGS
# ============================================================

# Estimate FPGS using random forest

fit_cont_rf <- fpgs(
  data = df_cont,
  outcome = "lead",
  treatment = "smoking",
  covariates = covariates,
  outcome_type = "continuous",
  method = "random_forest",
  folds = 5
)

fit_cont_rf



# Checking for random-forest FPGS components


summary(fit_cont_rf$mu0_hat)

summary(fit_cont_rf$mu1_hat)

length(fit_cont_rf$mu0_hat)

length(fit_cont_rf$mu1_hat)

anyNA(fit_cont_rf$mu0_hat)

anyNA(fit_cont_rf$mu1_hat)


# 17. ATE estimates using random-forest FPGS

ate_cont_rf <- fpgs_ate(
  fit = fit_cont_rf,
  variance = "none"
)

ate_cont_rf


# COMPARING LINEAR REGRESSION AND RANDOM FOREST

# Comparison table


comparison_cont <- data.frame(
  estimator = ate_cont$estimator,
  linear_regression = ate_cont$estimate,
  random_forest = ate_cont_rf$estimate
)

comparison_cont

# COMPUTATION TIME FOR BOOTSTRAP INFERENCE

B_time <- 1000

# Bootstrap inference: linear-regression FPGS


time_linear_boot <- system.time({

  ate_cont_boot_time <- fpgs_ate(
    fit = fit_cont,
    variance = "bootstrap",
    B = B_time,
    conf.level = 0.95
  )

})

# Results and computation time


ate_cont_boot_time

time_linear_boot

cat(
  "\nElapsed time for linear-regression bootstrap:",
  time_linear_boot["elapsed"],
  "seconds\n"
)


# ============================================================
# RANDOM-FOREST FPGS
# ============================================================


# ------------------------------------------------------------
# Bootstrap inference: random-forest FPGS
# ------------------------------------------------------------

time_rf_cont_boot <- system.time({

  ate_cont_rf_boot_time <- fpgs_ate(
    fit = fit_cont_rf,
    variance = "bootstrap",
    B = B_time,
    conf.level = 0.95
  )

})


# ------------------------------------------------------------
# Results and computation time
# ------------------------------------------------------------

ate_cont_rf_boot_time

time_rf_cont_boot

cat(
  "\nElapsed time for random-forest bootstrap:",
  time_rf_cont_boot["elapsed"],
  "seconds\n"
)


#Timing comparison

timing_comparison_cont <- data.frame(
  outcome = "Continuous (blood lead)",
  method = c(
    "Linear regression",
    "Random forest"
  ),
  B = B_time,
  elapsed_seconds = c(
    unname(time_linear_boot["elapsed"]),
    unname(time_rf_cont_boot["elapsed"])
  )
)

timing_comparison_cont$elapsed_minutes <-
  timing_comparison_cont$elapsed_seconds / 60

timing_comparison_cont

# Relative computation time


rf_time_ratio_cont <-
  unname(time_rf_cont_boot["elapsed"]) /
  unname(time_linear_boot["elapsed"])

cat(
  "\nRandom-forest bootstrap took",
  round(rf_time_ratio_cont, 2),
  "times as long as the linear-regression bootstrap.\n"
)


