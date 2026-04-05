
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
    month = factor(month, levels = 1:12, ordered = TRUE),
    dow = factor(dow, levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"),
                 ordered = TRUE),
    trend = as.integer(date - as.Date("2020-01-01")),
    #is_covid = ifelse(year %in% 2020:2021, "COVID Period", "Post-COVID")
    
  )

cycle_cv_df <- cycle_daily_df %>%
  mutate(month = as.numeric(month))

# 3. Model Fitting
# Note: Use factor(month) in formulas for M1-M3
model_formulas <- list(
  M0 = count ~ temp_mean + weekend + month_num,
  M1 = count ~ temp_mean + weekend + trend + factor(month) + factor(dow),
  M2 = count ~ weekend + trend + factor(month) + factor(dow) +
    temp_mean + I(temp_mean^2),
  M3 = log(count + 1) ~ weekend + trend + factor(month) + factor(dow) +
    temp_mean + I(temp_mean^2) + is_covid
)

m0 <- lm(model_formulas[["M0"]], data = cycle_daily_df)
m1 <- lm(model_formulas[["M1"]], data = cycle_daily_df)
m2 <- lm(model_formulas[["M2"]], data = cycle_daily_df)
m3 <- lm(model_formulas[["M3"]], data = cycle_daily_df)

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
    
    # Case for M3 
    if (model_name == "M3") {
      mu_log    <- mu
      sigma_log <- sigma
      
      # Back-transform to count scale
      mu_count    <- exp(mu_log) - 1
      sigma_count <- exp(mu_log) * sigma_log  # delta method
      #sigma_count <- sqrt((exp(sigma_log^2) - 1) * exp(2 * mu_log + sigma_log^2))
      
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
    
    mean_M3_count    = exp(mean_M3) - 1,
    sigma_count = exp(mean_M3) * sd_M3,
    
    se_M3  = (count - mean_M3_count)^2,
    ds_M3  = (count - mean_M3_count)^2 / sigma_count^2 + 2 * log(sigma_count),
    
    lwr_M3 = exp(mean_M3 - qnorm(0.975) * sd_M3) - 1,
    upr_M3 = exp(mean_M3 + qnorm(0.975) * sd_M3) - 1,
    
    is_M3  = (upr_M3 - lwr_M3) +
      (2 / 0.05) * pmax(lwr_M3 - count, 0) +
      (2 / 0.05) * pmax(count - upr_M3, 0)
    
    
    #mean_M3_count = exp(mean_M3) - 1,
    #se_M3 = (count - mean_M3_count)^2,
    #ds_M3 = (log(count + 1) - mean_M3)^2 / sd_M3^2 + 2 * log(sd_M3),
    #mu_count    = exp(mean_M3) - 1,
    #sigma_count = exp(mean_M3) * sd_M3,
    #ds_M3 = (count - mu_count)^2 / sigma_count^2 + 2 * log(sigma_count),
    
    #lwr_M3 = exp(mean_M3 - qnorm(0.975) * sd_M3) - 1,
    #upr_M3 = exp(mean_M3 + qnorm(0.975) * sd_M3) - 1,
    
    #lwr_M3 = mean_M3 - qnorm(0.975) * sd_M3,
    #upr_M3 = mean_M3 + qnorm(0.975) * sd_M3,
    #is_M3 = (upr_M3 - lwr_M3) +
      #(2 / 0.05) * pmax(lwr_M3 - log(count + 1), 0) +
      #(2 / 0.05) * pmax(log(count + 1) - upr_M3, 0)
  )

# Summarise scores by month
cv_month_allmodels <- score_month %>%
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

print(cv_month_allmodels)

#Model3 

cv_table2 <- score_month %>%
  group_by(month) %>%
  summarise(
    RMSE = round(sqrt(mean(se_M3)), 1),
    DS   = round(mean(ds_M3), 2),
    .groups = "drop"
  ) %>%
  mutate(
    Month = month.abb[as.integer(month)]
  ) %>%
  select(Month, RMSE, DS)

print(cv_table2)

# ...


# 7. Exploratory summaries and plots

summary_stats <- data.frame(
  Variable = c("Daily Count", "Mean Temp (°C)", "Min Temp (°C)", "Max Temp (°C)"),
  Mean = round(c(mean(cycle_daily_df$count),
                 mean(cycle_daily_df$temp_mean),
                 mean(cycle_daily_df$temp_min),
                 mean(cycle_daily_df$temp_max)), 1),
  SD = round(c(sd(cycle_daily_df$count),
               sd(cycle_daily_df$temp_mean),
               sd(cycle_daily_df$temp_min),
               sd(cycle_daily_df$temp_max)), 1),
  Min = round(c(min(cycle_daily_df$count),
                min(cycle_daily_df$temp_mean),
                min(cycle_daily_df$temp_min),
                min(cycle_daily_df$temp_max)), 1),
  Median = round(c(median(cycle_daily_df$count),
                   median(cycle_daily_df$temp_mean),
                   median(cycle_daily_df$temp_min),
                   median(cycle_daily_df$temp_max)), 1),
  Max = round(c(max(cycle_daily_df$count),
                max(cycle_daily_df$temp_mean),
                max(cycle_daily_df$temp_min),
                max(cycle_daily_df$temp_max)), 1)
)

print(summary_stats)

#Covid Indicator
cycle_daily_df <- cycle_daily_df %>%
  mutate(
    is_covid = ifelse(year %in% 2020:2021, "COVID Period", "Post-COVID")
  )


#EDA Plots
eda_time_series_plot <- ggplot(cycle_daily_df, aes(x = date, y = count)) +
  geom_line(alpha = 0.3, colour = "steelblue") +
  geom_smooth(colour = "darkred", se = FALSE) +
  labs(
    title = "Daily Cyclist Counts in Edinburgh, 2020-2025",
    x = "Date",
    y = "Daily Count"
  ) +
  theme_minimal()

eda_month_boxplot <- ggplot(cycle_daily_df, aes(x = month, y = count)) +
  geom_boxplot(fill = "blue", alpha = 0.6) +
  labs(
    title = "Count by Month",
    x = "Month",
    y = "Daily Count"
  ) +
  theme_minimal()

eda_dow_boxplot <- ggplot(cycle_daily_df, aes(x = dow, y = count)) +
  geom_boxplot(fill = "red", alpha = 0.6) +
  labs(
    title = "Count by Day of Week",
    x = "Day of Week",
    y = "Daily Count"
  ) +
  theme_minimal()

eda_temp_scatter <- ggplot(cycle_daily_df, aes(x = temp_mean, y = count)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess", colour = "red", se = FALSE) +
  labs(
    title = "Daily Cyclist Count vs Mean Temperature",
    x = "Mean Temperature (°C)",
    y = "Daily Count"
  ) +
  theme_minimal()

eda_covid_temp_plot <- ggplot(cycle_daily_df, aes(x = temp_mean, y = count, colour = is_covid)) +
  geom_point(alpha = 0.2) +
  geom_smooth(se = FALSE) +
  scale_colour_manual(values = c("COVID Period" = "darkred",
                                 "Post-COVID" = "steelblue")) +
  labs(
    title = "Cyclist Count vs Temperature by COVID Period",
    x = "Mean Temperature (°C)",
    y = "Daily Count",
    colour = ""
  ) +
  theme_minimal()

print(eda_time_series_plot)
print(eda_month_boxplot)
print(eda_dow_boxplot)
print(eda_temp_scatter)
print(eda_covid_temp_plot)

# 8 Baseline model outputs

m0_coef_table <- broom::tidy(m0) %>%
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 2),
    p.value = round(p.value, 4)
  )

print(m0_coef_table)

m0_diag_data <- data.frame(
  fitted = fitted(m0),
  resid = resid(m0),
  std_resid = rstandard(m0)
)

print(m0_diag_data)

m0_resid_fitted_plot <- ggplot(m0_diag_data, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  labs(
    title = "M0: Residuals vs Fitted",
    x = "Fitted values",
    y = "Residuals"
  ) +
  theme_minimal()

print(m0_resid_fitted_plot)

m0_qq_plot <- ggplot(m0_diag_data, aes(sample = std_resid)) +
  stat_qq(alpha = 0.4) +
  stat_qq_line(colour = "red") +
  labs(
    title = "M0: Normal Q-Q Plot",
    x = "Theoretical Quantiles",
    y = "Standardised Residuals"
  ) +
  theme_minimal()

print(m0_qq_plot)

# 9. Final model outputs (M3)

m3_coef_table <- broom::tidy(m3) %>%
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 2),
    p.value = round(p.value, 4)
  )

print(m3_coef_table)

m1_diag_data <- data.frame(
  fitted = fitted(m1),
  resid = resid(m1),
  std_resid = rstandard(m1)
)

m2_diag_data <- data.frame(
  fitted = fitted(m2),
  resid = resid(m2),
  std_resid = rstandard(m2)
)

m3_diag_data <- data.frame(
  fitted = fitted(m3),
  resid = resid(m3),
  std_resid = rstandard(m3)
)

m1_resid_fitted_plot <- ggplot(m1_diag_data, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  labs(title = "M1: Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
  theme_minimal()

m1_qq_plot <- ggplot(m1_diag_data, aes(sample = std_resid)) +
  stat_qq(alpha = 0.4) +
  stat_qq_line(colour = "red") +
  labs(title = "M1: Normal Q-Q Plot",
       x = "Theoretical Quantiles", y = "Standardised Residuals") +
  theme_minimal()

m2_resid_fitted_plot <- ggplot(m2_diag_data, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  labs(title = "M2: Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
  theme_minimal()

m2_qq_plot <- ggplot(m2_diag_data, aes(sample = std_resid)) +
  stat_qq(alpha = 0.4) +
  stat_qq_line(colour = "red") +
  labs(title = "M2: Normal Q-Q Plot",
       x = "Theoretical Quantiles", y = "Standardised Residuals") +
  theme_minimal()

m3_resid_fitted_plot <- ggplot(m3_diag_data, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  labs(title = "M3: Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
  theme_minimal()

m3_qq_plot <- ggplot(m3_diag_data, aes(sample = std_resid)) +
  stat_qq(alpha = 0.4) +
  stat_qq_line(colour = "red") +
  labs(title = "M3: Normal Q-Q Plot",
       x = "Theoretical Quantiles", y = "Standardised Residuals") +
  theme_minimal()

print(m1_resid_fitted_plot)
print(m1_qq_plot)
print(m2_resid_fitted_plot)
print(m2_qq_plot)
print(m3_resid_fitted_plot)
print(m3_qq_plot)

#!! !!! !! !!! !! ADD REMAINING DIAGNOSTIC PLOTS FOR M1 M2 M3 ??? !! !!! !! !!!

# 10. Table 3: trend coefficients from M1-M3

m1_daily_change <- coef(m1)["trend"]
m2_daily_change <- coef(m2)["trend"]

start_data <- cycle_daily_df %>%
  filter(date == as.Date("2020-01-01")) %>%
  slice(1)

start_pred_M3 <- exp(predict(m3, newdata = start_data)) - 1
m3_daily_change <- (start_pred_M3 + 1) * coef(m3)["trend"]

trend_coefs <- data.frame(
  Model = c("M1", "M2", "M3"),
  `Daily change (cyclists/day)` = round(c(m1_daily_change, m2_daily_change, m3_daily_change), 3),
  `Annual change (cyclists/year)` = round(c(m1_daily_change, m2_daily_change, m3_daily_change) * 365, 0)
)

print(trend_coefs)

# 10b. Trend significance for M1-M3

trend_significance <- data.frame(
  Model = c("M1", "M2", "M3"),
  Estimate = c(
    summary(m1)$coefficients["trend", "Estimate"],
    summary(m2)$coefficients["trend", "Estimate"],
    summary(m3)$coefficients["trend", "Estimate"]
  ),
  Std_Error = c(
    summary(m1)$coefficients["trend", "Std. Error"],
    summary(m2)$coefficients["trend", "Std. Error"],
    summary(m3)$coefficients["trend", "Std. Error"]
  ),
  t_value = c(
    summary(m1)$coefficients["trend", "t value"],
    summary(m2)$coefficients["trend", "t value"],
    summary(m3)$coefficients["trend", "t value"]
  ),
  p_value = c(
    summary(m1)$coefficients["trend", "Pr(>|t|)"],
    summary(m2)$coefficients["trend", "Pr(>|t|)"],
    summary(m3)$coefficients["trend", "Pr(>|t|)"]
  )
) %>%
  mutate(
    Estimate = round(Estimate, 4),
    Std_Error = round(Std_Error, 4),
    t_value = round(t_value, 2),
    p_value = round(p_value, 4)
  )

print(trend_significance)

# 10c. Annual growth as percentage of predicted count on 2020-01-01

start_data <- cycle_daily_df %>%
  filter(date == as.Date("2020-01-01")) %>%
  slice(1)

start_pred_M1 <- predict(m1, newdata = start_data)
start_pred_M2 <- predict(m2, newdata = start_data)
start_pred_M3 <- exp(predict(m3, newdata = start_data)) - 1

m1_daily_change <- coef(m1)["trend"]
m2_daily_change <- coef(m2)["trend"]

m3_trend_log <- coef(m3)["trend"]
m3_daily_change <- (start_pred_M3 + 1) * m3_trend_log

growth_rate_table <- data.frame(
  Model = c("M1", "M2", "M3"),
  Predicted_start_count = c(start_pred_M1, start_pred_M2, start_pred_M3),
  Daily_change = c(m1_daily_change, m2_daily_change, m3_daily_change)
) %>%
  mutate(
    Annual_change = Daily_change * 365,
    Growth_percent = 100 * Annual_change / Predicted_start_count,
    Predicted_start_count = round(Predicted_start_count, 1),
    Daily_change = round(Daily_change, 3),
    Annual_change = round(Annual_change, 0),
    Growth_percent = round(Growth_percent, 2)
  ) %>%
  select(Model, Predicted_start_count, Daily_change, Annual_change, Growth_percent)

print(growth_rate_table)

# 11. Section 5.1: compare temperature variables in M2

temp_model_formulas <- list(
  temp_mean = count ~ weekend + trend + factor(month) + factor(dow) +
    temp_mean + I(temp_mean^2),
  temp_min = count ~ weekend + trend + factor(month) + factor(dow) +
    temp_min + I(temp_min^2),
  temp_max = count ~ weekend + trend + factor(month) + factor(dow) +
    temp_max + I(temp_max^2)
)

temp_cv_results <- bind_rows(lapply(names(temp_model_formulas), function(temp_name) {
  bind_rows(lapply(years, function(yr) {
    
    test_df  <- cycle_cv_df %>% filter(year == yr)
    train_df <- cycle_cv_df %>% filter(year != yr)
    
    fit <- lm(temp_model_formulas[[temp_name]], data = train_df)
    pred <- predict(fit, newdata = test_df, se.fit = TRUE)
    
    mu <- pred$fit
    sigma <- sqrt(summary(fit)$sigma^2 + pred$se.fit^2)
    y <- test_df$count
    scores <- calc_scores(y, mu, sigma)
    
    data.frame(
      temp_var = temp_name,
      year = yr,
      RMSE = scores$RMSE,
      DS = mean(scores$DS),
      IS = mean(scores$IS)
    )
  }))
}))

temp_compare_table <- temp_cv_results %>%
  group_by(temp_var) %>%
  summarise(
    RMSE = round(mean(RMSE), 1),
    DS   = round(mean(DS), 2),
    IS   = round(mean(IS), 1),
    .groups = "drop"
  )



# 12. Marginal temperature effects for chosen temperature variable

print(temp_compare_table)

m2_temp_min <- lm(
  count ~ weekend + trend + factor(month) + factor(dow) +
    temp_min + I(temp_min^2),
  data = cycle_daily_df
)

beta1 <- coef(m2_temp_min)["temp_min"]
beta2 <- coef(m2_temp_min)["I(temp_min^2)"]

marginal_effect_5 <- beta1 + 2 * beta2 * 5
marginal_effect_15 <- beta1 + 2 * beta2 * 15

temp_effect_table <- data.frame(
  Temperature = c("5C", "15C"),
  `Marginal effect (cyclists per +1C)` = round(c(marginal_effect_5, marginal_effect_15), 1)
)

print(temp_effect_table)

# 13. COVID indicator extension check

cycle_daily_df <- cycle_daily_df %>%
  mutate(is_covid_binary = ifelse(year %in% 2020:2021, 1, 0))

m3_covid <- lm(
  log(count + 1) ~ weekend + trend + factor(month) + factor(dow) +
    temp_mean + I(temp_mean^2) + is_covid_binary,
  data = cycle_daily_df
)

start_pred_M3 <- exp(predict(m3, newdata = start_data)) - 1

start_data <- cycle_daily_df %>% filter(date == as.Date("2020-01-01")) %>% slice(1)

start_pred_M3_covid <- exp(predict(m3_covid, newdata = start_data)) - 1

m3_daily_change <- (start_pred_M3 + 1) * coef(m3)["trend"]
m3_covid_daily_change <- (start_pred_M3_covid + 1) * coef(m3_covid)["trend"]

trend_comparison <- data.frame(
  Model = c("M3", "M3 + COVID"),
  Daily_change = round(c(m3_daily_change, m3_covid_daily_change), 3),
  Annual_change = round(c(m3_daily_change, m3_covid_daily_change) * 365, 0)
)

print(trend_comparison)

# 14. Extrapolation to 2035 under fixed summer weekday conditions

future_dates <- data.frame(
  date = seq(as.Date("2020-01-01"), as.Date("2035-12-31"), by = "month")
) %>%
  mutate(
    year = year(date),
    month = factor(7, levels = 1:12, ordered = TRUE),   # July
    dow = factor("Wed", levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"),
                 ordered = TRUE),
    weekend = 0,
    temp_mean = mean(cycle_daily_df$temp_mean, na.rm = TRUE),
    trend = as.integer(date - as.Date("2020-01-01")),
    #is_covid = ifelse(year %in% 2020:2021, "COVID Period", "Post-COVID")
  )

future_dates$pred_log <- predict(m3, newdata = future_dates)
future_dates$pred_count <- exp(future_dates$pred_log) - 1

extrapolation_plot <- ggplot(future_dates, aes(x = date, y = pred_count)) +
  geom_line(colour = "darkblue") +
  geom_hline(yintercept = 15000, linetype = "dashed", colour = "red") +
  labs(
    title = "Projected Cycling Demand to 2035 under Summer Weekday Conditions",
    x = "Date",
    y = "Predicted Daily Count"
  ) +
  theme_minimal()

print(extrapolation_plot)

threshold_15000 <- future_dates %>%
  filter(pred_count >= 15000) %>%
  slice(1) %>%
  select(date, pred_count)

print(threshold_15000)

# 15. Poorly fit period for final model

cycle_daily_df$m3_fitted_log <- fitted(m3)
cycle_daily_df$m3_fitted_count <- exp(cycle_daily_df$m3_fitted_log) - 1
cycle_daily_df$m3_resid_count <- cycle_daily_df$count - cycle_daily_df$m3_fitted_count

resid_time_plot <- ggplot(cycle_daily_df, aes(x = date, y = m3_resid_count, colour = factor(year))) +
  geom_line() +
  labs(
    title = "Final Model Residuals over Time",
    x = "Date",
    y = "Residual (count scale)",
    colour = "Year"
  ) +
  theme_minimal()

print(resid_time_plot)

poor_period_table <- cycle_daily_df %>%
  mutate(year_month = format(date, "%Y-%m")) %>%
  group_by(year_month) %>%
  summarise(
    mean_residual = mean(m3_resid_count),
    RMSE = sqrt(mean(m3_resid_count^2)),
    .groups = "drop"
  ) %>%
  arrange(desc(RMSE))

print(head(poor_period_table, 10))

# 16. Original plot: monthly prediction bias over time

bias_plot_data <- cycle_daily_df %>%
  mutate(year_month = format(date, "%Y-%m")) %>%
  group_by(year_month) %>%
  summarise(
    mean_residual = mean(m3_resid_count),
    .groups = "drop"
  )

original_bias_plot <- ggplot(bias_plot_data, aes(x = year_month, y = mean_residual, group = 1)) +
  geom_line(colour = "purple") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  labs(
    title = "Monthly Mean Prediction Bias of Final Model",
    x = "Year-Month",
    y = "Mean Residual"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

print(original_bias_plot)

