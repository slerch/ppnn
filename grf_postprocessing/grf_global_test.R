rm(list=ls())

library(grf)
?quantile_forest

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_all.Rdata"))

## ignore sm_mean and sm_var due to missing and NA values leading to problems in estimation
data$sm_mean <- NULL
data$sm_var <- NULL

library(quantregForest)
library(scoringRules)
library(lubridate)

# train and eval dates (any date in date = valid date of fc and!)
start_eval <- as.Date("2016-01-01 00:00", tz = "UTC")
end_eval <- as.Date("2016-12-31 00:00", tz = "UTC") 

# prepare evaluation data set for faster subsetting later on
data_eval_all <- subset(data, date >= start_eval & date <= end_eval)

end_train <- start_eval - days(2)
start_train <- as.Date("2015-01-01 00:00", tz  = "UTC")

# prepare training data set collection for faster subsetting later on
data_train_all <- subset(data, date >= start_train & date <= end_train)

qt_levels <- seq(1/51, 50/51, by = 1/51)
qts_save <- matrix(NA, nrow = nrow(data_eval_all), ncol = length(qt_levels))

y_train <- data_train_all$obs
X_train <- data_train_all[!(names(data_train_all) %in% c("obs", "date", "station"))]

# estimate model on training data
grF_model <- quantile_forest(X_train, 
                             y_train,
                             quantiles = qt_levels,
                             num.trees = 500,
                             num.threads = 4)

X_eval <- data_eval_all[!(names(data_eval_all) %in% c("obs", "date", "station"))]

grF_prediction <-   predict(grF_model,
                            X_eval,
                            quantiles = qt_levels)


