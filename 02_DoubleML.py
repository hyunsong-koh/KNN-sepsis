import os
import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestRegressor
import doubleml as dml
from doubleml import DoubleMLData

# Reproducibility details
# Python version 3.9.23 was used. “doubleml”, version 0.10.1. was used.

# Define working directory
os.chdir("YOUR FOLDER")

# import csv file
df = pd.read_csv("data/knn_data_toy.1k.csv")

# # Define outcome, treatment variables and covariates
Y_COL = "seps"
T_COL = "invfpod_log"
X_COLS = [c for c in df.columns if c not in [Y_COL, T_COL]]

# 1. Set DoublML DATA
dml_data = DoubleMLData(
   df, y_col=Y_COL, d_cols=T_COL, x_cols=X_COLS
)

# 2. ML model

# Build ML model
ml_model = RandomForestRegressor(
    n_estimators=500, max_features=20, max_depth=5, min_samples_leaf=2, random_state=42
)

# 3. APE estimation

# Set random seed for reproducibility
np.random.seed(42)

# Set up DoubleML model
dml_model = dml.DoubleMLPLR(dml_data, ml_model, ml_model, n_rep=3)

# Model Fitting
print(dml_model.fit())

# Resutls of Average Partial Effect (APE) estimation
summary = dml_model.summary
print("Model Summary:\n", summary)

# 4. Sensitivity Analysis
dml_model.sensitivity_analysis()
print(dml_model.sensitivity_summary)
