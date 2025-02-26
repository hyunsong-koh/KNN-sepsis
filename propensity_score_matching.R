library(tidyverse)    
library(MatchIt)        
library(cobalt)         
library(rbounds)        
library(marginaleffects)

# Reproducibility details
# R version 4.4.1 was used.

# Set working directory (replace with your actual folder path)
setwd("YOUR FOLDER")

# Import CSV data
data <- read_csv("Your_DATA")

# Create Treatment variable
median_invfpod_rt <- median(data$invfpod_rt, na.rm=TRUE)
data$treatment <- ifelse(data$invfpod_rt > median_invfpod_rt, 1, 0)

# Remove 'invfpod_rt'
data <- data %>% dplyr::select(-invfpod_rt)

# Perform propensity score matching (exclude the outcome variable 'seps_v2' from matching)
set.seed(0)
m1 <- matchit(treatment ~ . -seps_v2, 
              data=data, 
              method="nearest", # Nearest neighbor matching
              distance="glm", 
              ratio=1, 
              replace=FALSE, 
              caliper=0.2)


# Check standardized mean differences (SMD) before and after matching
bal.tab(m1, un=TRUE)

# Visualize SMD using love.plot to assess balance after matching
love.plot(m1, thresholds=0.1,
          drop.distance=TRUE, 
          binary='std',
          position='bottom',
          var.order='unadjusted',
          size=2
)

# Create matched dataset
matched_df <- match.data(m1)

# Estimate ATT
lr <- glm(seps_v2 ~ . - PS - distance - weights - subclass, # remove the matching-related variables (PS, distance, weights, subclass)
                  family = binomial, # using logistic regression
                  data=matched_df)

lr_compare = avg_comparisons(lr, 
                             variables="treatment",
                             vcov=~subclass,
                             newdata = subset(matched_df, treatment==1))

# ATT estimate results
print(lr_compare) 

# ATT sensitivity test (psens only 1:1)
x <- matched_df %>% filter(treatment == 1) %>% pull(seps_v2)  
y <- matched_df %>% filter(treatment == 0) %>% pull(seps_v2)  
psens_sensitivity_analysis1 <- psens(x = x, y = y, Gamma = 3, GammaInc = 0.1)  # Adjust Gamma and GammaInc values as needed for your analysis

# Summary of sensitivity analysis results
summary(psens_sensitivity_analysis1)