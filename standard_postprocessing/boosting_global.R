## implementation of EMOS model with boosting 
## training on data from 2015 only
## here: results for some selected tuning parameters
## here: based on prepared Rdata file

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
train_start <- as.Date("2015-01-01 00:00 UTC") 

data_train <- subset(data, date >= train_start & date <= train_end)
data_train <- data_train[complete.cases(data_train), ]

crch_model <- crch(obs ~ .|.,
                   data = data_train[,-which(names(data) %in% c("date", "station"))],
                   dist = "gaussian",
                   link.scale = "log",
                   method = "boosting",
                   mstop = "aic")


eval_start <- as.Date("2016-01-01 00:00 UTC")
eval_end <- as.Date("2016-12-31 00:00 UTC")
eval_dates <- seq(eval_start, eval_end, by = "1 day")

data_eval_all <- subset(data, date >= eval_start & date <= eval_end)

loc <- as.numeric(predict(crch_model, data_eval_all, type = "location"))
sc <- as.numeric(predict(crch_model, data_eval_all, type = "scale"))

crps_boost <- crps_norm(y = data_eval_all$obs, mean = loc, sd = sc)
summary(crps_boost)

## tuning parameters for boosting:
# maxit = number of boosting iterations, default = 100
# nu = boosting step size
# mstop = c("max", "aic", "bic", "cv")
# nfolds = number of folds if mstop = cv
# link.scale = c("log", "identity", "quadratic")

# test in cluster experiments