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
data_train_all <- data_train_all[complete.cases(data_train_all),] 

nQ <- 10 # number of quantiles
qt_levels <- seq(1/(nQ+1), nQ/(nQ+1), by = 1/(nQ+1))
qts_save <- matrix(NA, nrow = nrow(data_eval_all), ncol = length(qt_levels))

stations <- unique(data$station)

for(this_station in stations){
  # progress indicator
  progind <- which(stations == this_station)
  if(progind %% 10 == 0){
    cat(progind, "of", length(stations), "started at", paste(Sys.time()), "\n")
  }
  
  data_train_thisstation <- subset(data_train_all, station == this_station)
  data_train <- data_train_thisstation[complete.cases(data_train_thisstation),]
  y_train <- data_train$obs
  X_train <- data_train[!(names(data_train) %in% c("obs", "date", "station"))]
  
  # omit stations with too few observation to successfully estimate model
  if(nrow(data_train) <= 10){next}
  
  # estimate model on training data
  grF_model <- quantile_forest(X_train, 
                               y_train,
                               quantiles = qt_levels,
                               num.threads = 2)
  
  # compute quantiles on evaluation data
  data_eval_thisstation <- subset(data_eval_all, station == this_station)
  X_eval <- data_eval_thisstation[!(names(data_eval_thisstation) %in% c("obs", "date", "station"))]

  grF_prediction <-   predict(grF_model,
                              X_eval,
                              quantiles = qt_levels)
  
  ind_this_station_eval <- which(data_eval_all$station == this_station)
  
  qts_save[ind_this_station_eval,] <- grF_prediction
}

ind_noNAinRow <- which(!apply(qts_save, 1, anyNA))
ind_use <- intersect(which(!is.na(data_eval_all$obs)), ind_noNAinRow)

dim(qts_save)[1] - length(ind_use)

grf_crps <- crps_sample(y = data_eval_all$obs[ind_use],
                        dat = qts_save[ind_use,])

summary(grf_crps) # mean 1.0942 (with 1/51 quantiles)

## mean 0.9979 with qt_levels <- seq(1/11, 10/11, by = 1/11)
## mean 1.0177 with 1/21 quantile steps
## mean 1.0143 with 1/21 quantile steps, but default number of trees (2000)
## mean 0.9973 with 1/6 quantile steps
# 0.9945 with 10 quantiles and default option for trees (2000)

## further setting sample.fraction = 0.632, leads to a mean CRPS of 0.983 
## (suggested by Maxime Taillardat)