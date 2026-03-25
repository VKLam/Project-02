
# ============================================================
# Project 2: Modelling daily cycling demand in Edinburgh
# Authors: [Karthik Sutheeshkumar], [Kiet Lam], [Leon Godtfredsen]
# ============================================================

#
# Place the code needed in the Report_project02.Rmd, including documentation.
#

# 0. Packages
library(ggplot2); library(dplyr); library(tidyr)
library(lubridate); library(knitr); library(kableExtra)
library(patchwork); library(broom)

# 1. Load data
load('cycle_daily_df.Rdata')

# 2. Data Wrangling 
cycle_daily_df <- cycle_daily_df %>%
  mutate(
    month = factor(month, levels = 1:12, ordered = TRUE),
    
    # Task 2: dow with explicit levels 
    dow = factor(dow, levels = c("Sun","Mon","Tue","Wed","Thu","Fri","Sat"), ordered = TRUE),
    
    
    # Task 3: trend as integer days since 2020-01-01
    trend = as.integer(date - as.Date("2020-01-01"))
    # Task 1: month as ordered factor for plots/inference

    # Task 2: dow with explicit levels 

    # Task 3: trend as integer days since 2020-01-01
    
  )

# 3. Model Fitting
# Note: Use factor(month) in formulas for M1-M3
m0 <- lm(count ~ temp_mean + weekend + month, data = cycle_daily_df)
m1<-lm(count~ temp_mean + weekend +trend+factor(month)+factor(dow),data=cycle_daily_df)
m2<- lm(count~ I(temp_mean^2)+weekend +trend+factor(month)+factor(dow),data=cycle_daily_df)
plot(m2)
m3<- lm(log(count+1)~ I(temp_mean^2)+weekend +trend+factor(month)+factor(dow),data=cycle_daily_df)
# m1 <- ...
# m2 <- ...
# m3 <- ...

# 4. Cross-Validation Functions

calc_scores <- function(y, mu, sigma, alpha = 0.05) {
  # y     : vector of observed values
  # mu    : vector of predictive means (from predict(fit, newdata=test)$fit)
  # sigma : vector of predictive SDs. Combine residual error and mean uncertainty:
  #         sigma = sqrt(summary(fit)$sigma^2 + predict(fit, newdata=test, se.fit=TRUE)$se.fit^2)
  # Returns a named list with RMSE, MAE, DS, IS
  
  # Implement RMSE, MAE, DS, and IS here
  
  RMSE <- sqrt(mean((y-mu)^2))
  MAE <- mean(abs(y-mu))
  DS <- (y-mu)^2/sigma^2 + 2 * log(sigma)

}

# 5. Leave-One-Year-Out CV Loop 
# 6. CV by Month 
# ...
