## Implementation of quantile regression forest approach to post-processing
## proposed by Taillardat et al (2016, Monthly Weather Review)
## partly based on R code provided by Maxime Taillardat

rm(list=ls())

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_all.Rdata"))

## ignore sm_mean and sm_var due to missing and NA values leading to problems in estimation
data$sm_mean <- NULL
data$sm_var <- NULL

library(quantregForest)
library(scoringRules)
library(lubridate)

start_eval <- as.Date("2016-01-01 00:00", tz = "UTC")
end_eval <- as.Date("2016-12-31 00:00", tz = "UTC") 

data_eval_all <- subset(data, date >= start_eval & date <= end_eval)

end_train <- start_eval - hours(48)
start_train <- as.Date("2015-01-01 00:00", tz = "UTC")

qt_levels <- seq(1/51, 50/51, by = 1/51)
qts_save <- matrix(NA, 
                   nrow = nrow(data_eval_all), 
                   ncol = length(qt_levels))

data_train <- subset(data, date >= start_train & date <= end_train)
data_train <- data_train[complete.cases(data_train),]

data_eval <- subset(data, date >= start_eval & date <= end_eval)

# model
qrF_model <- quantregForest(x = data_train[, !(names(data_train) %in% c("date", "station", "obs"))], 
                            y = data_train$obs)

# compute quantiles on evaluation data
qrF_prediction <-   predict(qrF_model,
                            newdata = data_eval[, !(names(data_train) %in% c("date", "station", "obs"))],
                            what = qt_levels,
                            all = TRUE)

# qt_levels <- seq(1/51, 50/51, by = 1/51)
# qts_save <- matrix(NA, nrow = nrow(data_eval), ncol = length(qt_levels))
# 
# qrF_prediction <-   predict(qrF_model,
#                             data_eval,
#                             what = qt_levels,
#                             all = TRUE)

# Evaluation
ind_noNAinRow <- which(!apply(qrF_prediction, 1, anyNA))
ind_use <- intersect(which(!is.na(data_eval$obs)), ind_noNAinRow)

qrf_crps <- crps_sample(y = data_eval$obs[ind_use],
                        dat = qrF_prediction[ind_use,])

summary(qrf_crps)

# Variable importance
qrF_model$importance
importance(qrF_model, type=1, scale=TRUE)

plot(importance(qrF_model, type=1, scale=TRUE))
