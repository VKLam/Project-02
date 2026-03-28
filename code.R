
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
model_formulas <- list(
  M0 = count ~ temp_mean + weekend + factor(month),
  M1 = count ~ temp_mean + weekend + trend + factor(month) + factor(dow),
  M2 = count ~ temp_mean + I(temp_mean^2)  + weekend + trend + factor(month) + factor(dow),
  M3 = log(count+1) ~ temp_mean + I(temp_mean^2)+ weekend +trend +factor(month) + factor(dow)
)

m0 <- lm(model_formulas[["M0"]], data = cycle_daily_df)
m1<-lm(model_formulas[["M1"]], data=cycle_daily_df)
m2<- lm(model_formulas[["M2"]], data=cycle_daily_df)
m3<- lm(model_formulas[["M3"]], data=cycle_daily_df)

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
      y_log <- log(test_df$count + 1)  # log scale for DS and IS
      y     <- test_df$count            # count scale for RMSE and MAE
      
      mu_log <- mu                      # log scale
      mu_count <- exp(mu) - 1           # count scale
      
      sigma_log <- sigma                # log scale
      
      # RMSE and MAE on count scale
      RMSE <- sqrt(mean((y - mu_count)^2))
      MAE  <- mean(abs(y - mu_count))
      
      # DS and IS on log scale
      scores_log <- calc_scores(y_log, mu_log, sigma_log)
      
      return(data.frame(
        model = model_name,
        year  = yr,
        RMSE  = RMSE,
        MAE   = MAE,
        DS    = mean(scores_log$DS),
        IS    = mean(scores_log$IS)
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

# ...
