rm(list=ls())

library(grf)

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_all.Rdata"))
load(paste0(data_dir, "data_all_stationInfo.Rdata"))

data$station_alt <- station_alt
data$station_orog <- station_orog
data$station_lsm <- station_lsm
rm(station_alt, station_orog, station_lsm)

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
# data_train_all <- data_train_all[complete.cases(data_train_all),] 

# qt_levels <- seq(1/51, 50/51, by = 1/51)
nQ <- 20 # number of quantiles
qt_levels <- seq(1/(nQ+1), nQ/(nQ+1), by = 1/(nQ+1))

y_train <- data_train_all$obs
X_train <- data_train_all[!(names(data_train_all) %in% c("obs", "date", "station"))]

# estimate model on training data
grF_model <- quantile_forest(X_train, 
                             y_train,
                             quantiles = qt_levels,
                             num.threads = 2)
# 
# X_eval <- data_eval_all[!(names(data_eval_all) %in% c("obs", "date", "station"))]
# 
# grF_prediction <-   predict(grF_model,
#                             newdata = X_eval,
#                             quantiles = qt_levels)
# 
# ind_noNAinRow <- which(!apply(grF_prediction, 1, anyNA))
# ind_use <- intersect(which(!is.na(data_eval_all$obs)), ind_noNAinRow)
# 
# dim(grF_prediction)[1] - length(ind_use)
# 
# grf_crps <- crps_sample(y = data_eval_all$obs[ind_use],
#                         dat = grF_prediction[ind_use,])
# 
# summary(grf_crps)

##
## ?????
##
# 
# Fehler in quantile_train(quantiles, regression.splitting, data$default,  : 
#                            Not compatible with requested type: [type=character; target=double].