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
outcome_var <- "seps" # Binary outcome (sepsis)
treatment_var <- "invfpod_log" # Continuous treatment (parenteral nutrition duration)

# Create Formula for analysis
covariates <- setdiff(names(data), c(outcome_var, treatment_var))
outcome_formula <- as.formula(paste(outcome_var, "~", treatment_var))

# Method 1: IPW (Inverse Probability Weighting)
# 1.1 Calculate Weights
set.seed(123)
ipw_model <- ipwpoint(
  exposure = invfpod_log,
  family = "gaussian",
  numerator = ~1,
  denominator = ~ var1 + var2 + var3 + var4 + var5 + var6 + var7 + var8 + var9 + 
    var10 + var11 + var12 + var13 + var14 + var15 + var16 + var17 + 
    var18 + var19 + var20 + var21 + var22 + var23 + var24 + var25 + 
    var26 + var27 + var28 + var29,
  data = as.data.frame(data),
  trunc = 0.01 # Truncate extreme weights (1% top/bottom)
)

# 1.2 Add ipw_weights to data
data$ipw_weights <- ipw_model$weights.trunc

# 1.3 Love Plot: Checking Covariate Balance
plot_formula <- as.formula(paste("invfpod_log ~", paste(covariates, collapse = " + ")))

love.plot(
  x = plot_formula,
  data = data,
  weights = data$ipw_weights,
  stats = "cor",       # Correlation for continuous treatment
  abs = TRUE,
  un = TRUE,
  thresholds = c(cor = 0.2),
  var.order = "unadjusted",
  title = "Covariate Balance Check (Correlations)"
)

# 1.4 Regression Analysis (Weighted GLM)
des_ipw <- svydesign(ids = ~1, weights = ~ipw_weights, data = data)
fit_ipw <- svyglm(outcome_formula, design = des_ipw, family = gaussian(link = "identity"))
summary(fit_ipw)

# Method 2: WeightIt (CBPS)
# Delete ipw_weights
data <- subset(data, select=-c(ipw_weights))

# 2.1 Calculate Weights
set.seed(123)
w_cbps <- weightit(
  formula  = invfpod_log ~ .,
  data     = data[, c(treatment_var, covariates), drop = FALSE],
  method = "cbps",
  estimand = "ATE",
  trim     = c(0.01, 0.99)
)

# 2.2 Add weights to data
data$cbps_weights <- w_cbps$weights

# 2.3 Love Plot: Checking Covariate Balance
plot_formula <- as.formula(paste("invfpod_log ~", paste(covariates, collapse = " + ")))

love.plot(
  x = plot_formula,
  data = data,
  weights = data$cbps_weights,
  stats = "cor",       # Correlation for continuous treatment
  abs = TRUE,
  un = TRUE,
  thresholds = c(cor = 0.2),
  var.order = "unadjusted",
  title = "Covariate Balance Check (Correlations)"
)


# 2.4 Regrssion Analysis (Weighted GLM)
des_cbps <- svydesign(ids = ~1, weights = ~cbps_weights, data = data)
fit_cbps <- svyglm(outcome_formula, design = des_cbps, family = gaussian(link = "identity"))
summary(fit_cbps)

