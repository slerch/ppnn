## implementation of EMOS model with boosting 
## here: results for some selected tuning parameters
## here: based on prepared Rdata file
## "_alldata" indicates use of all available training data before 2016

## --------------- data preparation --------------- ##

rm(list=ls())

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_all.Rdata"))

# remove sm data (missing values)
data$sm_mean <- NULL
data$sm_var <- NULL
head(data)

library(scoringRules)
library(lubridate)
library(crch)

# determine training set
train_end <- as.Date("2016-01-01 00:00 UTC") - days(2)
train_start <- data$date[1]

data_train <- subset(data, date >= train_start & date <= train_end)
data_train <- data_train[complete.cases(data_train), ]

crch_model <- crch(obs ~ .|.,
                   data = data_train[,-which(names(data) %in% c("date", "station"))],
                   dist = "gaussian",
                   link.scale = "log",
                   maxit = 1000,
                   nu = 0.05,
                   mstop = "max",
                   method = "boosting")


eval_start <- as.Date("2016-01-01 00:00 UTC")
eval_end <- as.Date("2016-12-31 00:00 UTC")
eval_dates <- seq(eval_start, eval_end, by = "1 day")

data_eval_all <- subset(data, date >= eval_start & date <= eval_end)

loc <- as.numeric(predict(crch_model, data_eval_all, type = "location"))
sc <- as.numeric(predict(crch_model, data_eval_all, type = "scale"))

crps_boost <- crps_norm(y = data_eval_all$obs, mean = loc, sd = sc)
summary(crps_boost)
