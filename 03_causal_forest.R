library(grf)
library(ggplot2)
library(dplyr)
library(stats)

# Reproducibility details
# R version 4.4.1 was used. “grf”, version 2.3.2. was used.

# Define working directory
setwd("YOUR FOLDER")

# import csv file
data <- read.csv("data/knn_data_toy.1k.csv")

# Define outcome and treatment variables
outcome <- as.vector(data$seps) # binary outcome (sepsis)
treatment <- as.vector(data$invfpod_log) # continuous treatment (parenteral nutrition duration)

## 1. Build causal forest model (GRF)
# Define covariates
X <- data[, !(names(data) %in% c("seps", "invfpod_log"))]
Y <- outcome
W <- treatment

# Generate causal forest model
set.seed(20250101, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
cf <- causal_forest(X = as.matrix(X), Y = outcome, W = treatment, Y.hat = NULL, W.hat = NULL, tune.parameters = "all", num.trees = 2000, tune.num.trees = 100)

# Perform calibration test
cf.fit <- test_calibration(cf)
print(cf.fit)

# Calculate fit index
coef_values <- coef(cf.fit)
beta_ATE <- coef_values['mean.forest.prediction']
beta_ITE <- coef_values['differential.forest.prediction']
fit_index <- abs(1-beta_ATE) + abs(1-beta_ITE)
print(fit_index)


## 2. APE and IPE estimation
# Average Partial Effect (APE) estimation
ate.result <- average_treatment_effect(cf, method = "AIPW", num.trees.for.weights = 1000, target.sample = "overlap")
ate.est <- ate.result[1]
ate.se <- ate.result[2]
ate.tstat <- ate.est / ate.se
ate.pvalue <- 1.96 * (1 - pnorm(abs(ate.tstat)))
ate.summary <- c(estimate = ate.est, std.error = ate.se, t.stat = ate.tstat, p.value = ate.pvalue)
print(ate.summary)

# Individual Partial Effect (IPE) estimation
ite <- predict(cf)$predictions
data[, "ite"] <- ite
hist(data$ite)
