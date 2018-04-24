## Implementation of quantile regression forest approach to post-processing
## proposed by Taillardat et al (2016, Monthly Weather Review)
## partly based on R code provided by Maxime Taillardat
## here: local = one model for every station

rm(list=ls())

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

data_eval_all <- subset(data, date >= start_eval & date <= end_eval)

end_train <- start_eval - hours(48)
start_train <- data$date[1]

# quantile levels
qt_levels <- seq(1/51, 50/51, by = 1/51)
qts_save <- matrix(NA, 
                   nrow = nrow(data_eval_all), 
                   ncol = length(qt_levels))

stations <- unique(data$station)

for(this_station in stations){
  # progress indicator
  progind <- which(stations == this_station)
  if(progind %% 10 == 0){
    cat(progind, "of", length(stations), "started at", paste(Sys.time()), "\n")
  }
  
  data_train <- subset(data,
                       date >= start_train & date <= end_train & station == this_station)
 
  data_train <- data_train[complete.cases(data_train),]
  
  data_eval <- subset(data,
                      date >= start_eval & date <= end_eval & station == this_station)
  # need to delete cases with NA obs in eval data to avoid error in qrf package
  # (even though not needed in predict function...)
  data_eval$obs <- NULL
  
  # omit stations with too few observation to successfully estimate model
  if(nrow(data_train) <= 10){next}
  
  # estimate model on training data
  qrF_model <- quantregForest(x = data_train[, !(names(data_train) %in% c("date", "station", "obs"))], 
                              y = data_train$obs,
                              nodesize = 5,
                              ntree = 500)
  
  # compute quantiles on evaluation data
  qrF_prediction <-   predict(qrF_model,
                              newdata = data_eval,
                              what = qt_levels,
                              all = TRUE)
  
  
  # position of fc cases for eval data
  # to save qrF prediction to row in qts_save
  qts_save[which(data_eval_all$station == this_station),] <- qrF_prediction
}


ind_noNAinRow <- which(!apply(qts_save, 1, anyNA))
ind_use <- intersect(which(!is.na(data_eval_all$obs)), ind_noNAinRow)

dim(qts_save)[1] - length(ind_use) # all NA cases in all eval data
dim(qts_save)[1] - length(ind_use) - sum(is.na(data_eval_all$obs)) # NA cases due to model estimation in eval data

qrf_crps <- crps_sample(y = data_eval_all$obs[ind_use],
                        dat = qts_save[ind_use,])

summary(qrf_crps)

# 0.8215