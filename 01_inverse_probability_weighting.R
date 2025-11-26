library(ipw)
library(WeightIt)   
library(survey)     
library(cobalt)

# Reproducibility details
# R version 4.4
# ipw version 1.2.1.1
# WeightIt   1.4.0
# survey   4.4.2
# cobalt   4.6.1


# Define working directory
setwd("YOUR FOLDER")

# import csv file
data <- read.csv("data/knn_data_toy.1k.csv")

# Define outcome and treatment variables
outcome_var <- "seps" # binary outcome (sepsis)
treatment_var <- "invfpod_log" # continuous treatment (parenteral nutrition duration)

# Create Formula for Propensity Score Models (Treatment ~ Covariates)
covariates <- setdiff(names(data), c(outcome_var, treatment_var))
ps_formula <- as.formula(paste(treatment_var, "~", paste(covariates, collapse = "+")))
outcome_formula <- as.formula(paste(outcome_var, "~", treatment_var))

# Method 1: IPW (Inverse Probability Weighting)
ipw_model <- ipwpoint(
  exposure = data[[treatment_var]],
  family = "gaussian",
  numerator = ~1,
  denominator = ps_formula,
  data = data,
  trunc = 0.01 # Truncate extreme weights (1% top/bottom)
)

# 1.2 Add weights to data
data$ipw_weights <- ipw_model$weights.trunc

# 1.3 Outcome Analysis (Weighted GLM)
des_ipw <- svydesign(ids = ~1, weights = ~ipw_weights, data = data)
fit_ipw <- svyglm(outcome_formula, design = des_ipw, family = gaussian(link = "identity"))

# Method 2: WeightIt (CBPS)
# 5.1 Calculate Weights
w_cbps <- weightit(
  ps_formula,
  data = data,
  method = "cbps",
  estimand = "ATE",
  over = FALSE
)

# 2.2 Add weights to data
data$cbps_weights <- w_cbps$weights

# 2.3 Outcome Analysis (Weighted GLM)
des_cbps <- svydesign(ids = ~1, weights = ~cbps_weights, data = data)
fit_cbps <- svyglm(outcome_formula, design = des_cbps, family = gaussian(link = "identity"))

# Love Plot: Checking Covariate Balance
love.plot(
  ps_formula,
  data = data,
  weights = list("Select IPW or weightit_CBPS"),    # IPW = data$ipw_weights, weightit_CBPS = data$cbps_weights
  stats = "cor",       # Correlation for continuous treatment
  abs = TRUE,
  thresholds = c(cor = 0.1),
  var.order = "unadjusted",
  title = "Covariate Balance Check (Correlations)"
)

