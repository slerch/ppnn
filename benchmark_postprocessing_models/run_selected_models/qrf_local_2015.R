## run selected benchmark models
## for documented code see model run code in ../standard_postprocessing/ and ../qrf_postprocessing/
## see note.md

## training 2015
## QRF local

## chosen tuning parameters:
##    ntree = 1000, nodesize = 10, mtry = 25, replace = TRUE

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
start_train <- as.Date("2015-01-01 00:00", tz  = "UTC")

# quantile levels
qt_levels <- seq(1/51, 50/51, by = 1/51)
qts_save <- matrix(NA, 
                   nrow = nrow(data_eval_all), 
                   ncol = length(qt_levels))

stations <- unique(data$station)

t1 <- Sys.time()

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
                              nodesize = 10,
                              ntree = 1000,
                              mtry = 25,
                              replace = TRUE)
  
  # compute quantiles on evaluation data
  qrF_prediction <-   predict(qrF_model,
                              newdata = data_eval,
                              what = qt_levels,
                              all = TRUE)
  
  
  # position of fc cases for eval data
  # to save qrF prediction to row in qts_save
  qts_save[which(data_eval_all$station == this_station),] <- qrF_prediction
}

t2 <- Sys.time()
t2-t1
# run time ~ 28 minutes

qrf_pred <- qts_save

# ind_noNAinRow <- which(!apply(qrf_pred, 1, anyNA))
# ind_use <- intersect(which(!is.na(data_eval_all$obs)), ind_noNAinRow)
# crps_qrf <- crps_sample(data_eval_all$obs[ind_use], qrf_pred[ind_use,])
# summary(crps_qrf)

df_out <- data.frame(date = data_eval_all$date,
                     station_id = data_eval_all$station,
                     qrf_pred)
str(df_out)

save(qrf_pred,
     file = "/home/sebastian/Projects/PP_NN/code/results/Rdata_files/qrf_local_2015.Rdata")
