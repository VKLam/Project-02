
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
    month_num = month,
    # Task 1: month as ordered factor for plots/inference
    month = factor(month, levels = 1:12, ordered = TRUE),
    
    # Task 2: dow with explicit levels 
    dow = factor(dow, levels = c("Sun","Mon","Tue","Wed","Thu","Fri","Sat"), ordered = TRUE),
    
    
    # Task 3: trend as integer days since 2020-01-01
    trend = as.integer(date - as.Date("2020-01-01"))

  )



# 3. Model Fitting
# Note: Use factor(month) in formulas for M1-M3
model_formulas <- list(
  M0 = count ~ temp_mean + weekend + month_num,
  M1 = count ~ temp_mean + weekend + trend + factor(month) + factor(dow),
  M2 = count ~ temp_mean + I(temp_mean^2)  + weekend + trend + factor(month) + factor(dow),
  M3 = log(count+1) ~ temp_mean + I(temp_mean^2)+ weekend +trend +factor(month) + factor(dow)
)

m0 <- lm(model_formulas[["M0"]], data = cycle_daily_df)
m1 <-lm(model_formulas[["M1"]], data=cycle_daily_df)
m2 <- lm(model_formulas[["M2"]], data=cycle_daily_df)
m3 <- lm(model_formulas[["M3"]], data=cycle_daily_df)

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
  l  <- mu - qnorm(1 - alpha/2) * sigma
  u  <- mu + qnorm(1 - alpha/2) * sigma
  IS <- (u-l) + (2/alpha)*pmax(l-y,0) + (2/alpha)*pmax(y-u,0)
  return(list(RMSE = RMSE, MAE = MAE, DS = DS, IS = IS))

}

# 5. Leave-One-Year-Out CV Loop 
cycle_cv_df <- cycle_daily_df %>% mutate(month = as.numeric(month))

years <-unique(cycle_cv_df$year)

cv_results <- bind_rows(lapply(names(model_formulas), function(model_name) {
  bind_rows(lapply(years, function(yr) {
    
    # Splitting into train and test
    test_df  <- cycle_cv_df %>% filter(year == yr)
    train_df <- cycle_cv_df %>% filter(year != yr)
    
    # Fitting model on training data 
    fit <- lm(model_formulas[[model_name]], data = train_df)
    
    # Predicting on test data
    pred  <- predict(fit, newdata = test_df, se.fit = TRUE)
    mu    <- pred$fit
    sigma <- sqrt(summary(fit)$sigma^2 + pred$se.fit^2)
    
    # Observed values from test_df
    y <- test_df$count
    if (model_name == "M3") {
      mu_log    <- mu
      sigma_log <- sigma
      
      # Back-transform to count scale
      mu_count    <- exp(mu_log) - 1
      sigma_count <- exp(mu_log) * sigma_log  # delta method
      
      y <- test_df$count
      
      CV_scores <- calc_scores(y, mu_count, sigma_count)
      
      return(data.frame(
        model = model_name,
        year  = yr,
        RMSE  = CV_scores$RMSE,
        MAE   = CV_scores$MAE,
        DS    = mean(CV_scores$DS),
        IS    = mean(CV_scores$IS)
      ))
    }
    
    # CV score predictions against observed values
    CV_scores <- calc_scores(y, mu, sigma)
    
    data.frame(
      model = model_name,
      year  = yr,
      RMSE  = CV_scores$RMSE,
      MAE   = CV_scores$MAE,
      DS    = mean(CV_scores$DS),
      IS    = mean(CV_scores$IS)
    )
  }))
}))

# Table 1: average scores across all held-out years
cv_table1 <- cv_results %>%
  group_by(model) %>%
  summarise(
    RMSE = round(mean(RMSE), 1),
    MAE  = round(mean(MAE), 1),
    DS   = round(mean(DS), 2),
    IS   = round(mean(IS),1),
    .groups = "drop"
  )

print(cv_table1)





# 6. CV by Month 


# Create an empty data frame to store observation-level CV predictions
cv_month_data <- cycle_cv_df %>%
  mutate(
    mean_M0 = NA_real_, sd_M0 = NA_real_,
    mean_M1 = NA_real_, sd_M1 = NA_real_,
    mean_M2 = NA_real_, sd_M2 = NA_real_,
    mean_M3 = NA_real_, sd_M3 = NA_real_
  )

# Fill in leave-one-year-out predictions for each model
for (yr in years) {
  
  train_data <- cycle_cv_df %>% filter(year != yr)
  test_ind   <- which(cycle_cv_df$year == yr)
  test_data  <- cycle_cv_df[test_ind, , drop = FALSE]
  
  # M0
  fit_M0 <- lm(model_formulas[["M0"]], data = train_data)
  pred_M0 <- predict(fit_M0, newdata = test_data, se.fit = TRUE)
  cv_month_data[test_ind, "mean_M0"] <- pred_M0$fit
  cv_month_data[test_ind, "sd_M0"] <- sqrt(pred_M0$se.fit^2 + summary(fit_M0)$sigma^2)
  
  # M1
  fit_M1 <- lm(model_formulas[["M1"]], data = train_data)
  pred_M1 <- predict(fit_M1, newdata = test_data, se.fit = TRUE)
  cv_month_data[test_ind, "mean_M1"] <- pred_M1$fit
  cv_month_data[test_ind, "sd_M1"] <- sqrt(pred_M1$se.fit^2 + summary(fit_M1)$sigma^2)
  
  # M2
  fit_M2 <- lm(model_formulas[["M2"]], data = train_data)
  pred_M2 <- predict(fit_M2, newdata = test_data, se.fit = TRUE)
  cv_month_data[test_ind, "mean_M2"] <- pred_M2$fit
  cv_month_data[test_ind, "sd_M2"] <- sqrt(pred_M2$se.fit^2 + summary(fit_M2)$sigma^2)
  
  # M3
  fit_M3 <- lm(model_formulas[["M3"]], data = train_data)
  pred_M3 <- predict(fit_M3, newdata = test_data, se.fit = TRUE)
  cv_month_data[test_ind, "mean_M3"] <- pred_M3$fit
  cv_month_data[test_ind, "sd_M3"] <- sqrt(pred_M3$se.fit^2 + summary(fit_M3)$sigma^2)
}

# Compute scores for each observation
score_month <- cv_month_data %>%
  mutate(
    se_M0 = (count - mean_M0)^2,
    ds_M0 = (count - mean_M0)^2 / sd_M0^2 + 2 * log(sd_M0),
    lwr_M0 = mean_M0 - qnorm(0.975) * sd_M0,
    upr_M0 = mean_M0 + qnorm(0.975) * sd_M0,
    is_M0 = (upr_M0 - lwr_M0) +
      (2 / 0.05) * pmax(lwr_M0 - count, 0) +
      (2 / 0.05) * pmax(count - upr_M0, 0),
    
    se_M1 = (count - mean_M1)^2,
    ds_M1 = (count - mean_M1)^2 / sd_M1^2 + 2 * log(sd_M1),
    lwr_M1 = mean_M1 - qnorm(0.975) * sd_M1,
    upr_M1 = mean_M1 + qnorm(0.975) * sd_M1,
    is_M1 = (upr_M1 - lwr_M1) +
      (2 / 0.05) * pmax(lwr_M1 - count, 0) +
      (2 / 0.05) * pmax(count - upr_M1, 0),
    
    se_M2 = (count - mean_M2)^2,
    ds_M2 = (count - mean_M2)^2 / sd_M2^2 + 2 * log(sd_M2),
    lwr_M2 = mean_M2 - qnorm(0.975) * sd_M2,
    upr_M2 = mean_M2 + qnorm(0.975) * sd_M2,
    is_M2 = (upr_M2 - lwr_M2) +
      (2 / 0.05) * pmax(lwr_M2 - count, 0) +
      (2 / 0.05) * pmax(count - upr_M2, 0),
    
    mean_M3_count = exp(mean_M3) - 1,
    se_M3 = (count - mean_M3_count)^2,
    ds_M3 = (log(count + 1) - mean_M3)^2 / sd_M3^2 + 2 * log(sd_M3),
    lwr_M3 = mean_M3 - qnorm(0.975) * sd_M3,
    upr_M3 = mean_M3 + qnorm(0.975) * sd_M3,
    is_M3 = (upr_M3 - lwr_M3) +
      (2 / 0.05) * pmax(lwr_M3 - log(count + 1), 0) +
      (2 / 0.05) * pmax(log(count + 1) - upr_M3, 0)
  )

# Summarise scores by month
cv_table2 <- score_month %>%
  group_by(month) %>%
  summarise(
    RMSE_M0 = round(sqrt(mean(se_M0)), 1),
    MAE_M0  = round(mean(abs(count - mean_M0)), 1),
    DS_M0   = round(mean(ds_M0), 2),
    IS_M0   = round(mean(is_M0), 1),
    
    RMSE_M1 = round(sqrt(mean(se_M1)), 1),
    MAE_M1  = round(mean(abs(count - mean_M1)), 1),
    DS_M1   = round(mean(ds_M1), 2),
    IS_M1   = round(mean(is_M1), 1),
    
    RMSE_M2 = round(sqrt(mean(se_M2)), 1),
    MAE_M2  = round(mean(abs(count - mean_M2)), 1),
    DS_M2   = round(mean(ds_M2), 2),
    IS_M2   = round(mean(is_M2), 1),
    
    RMSE_M3 = round(sqrt(mean(se_M3)), 1),
    MAE_M3  = round(mean(abs(count - mean_M3_count)), 1),
    DS_M3   = round(mean(ds_M3), 2),
    IS_M3   = round(mean(is_M3), 1)
  )

print(cv_table2)

# ...





# ...
