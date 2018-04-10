library(grf)
?quantile_forest


data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_prep_2015-16.Rdata"))

## ignore sm_mean and sm_var due to missing and NA values leading to problems in estimation
data_train$sm_mean <- NULL
data_train$sm_var <- NULL
data_eval$sm_mean <- NULL
data_eval$sm_var <- NULL

library(quantregForest)
library(scoringRules)

qt_levels <- seq(1/51, 50/51, by = 1/51)
qts_save <- matrix(NA, nrow = nrow(data_eval), ncol = length(qt_levels))

stations <- unique(data_eval$station)

for(this_station in stations){
  # progress indicator
  progind <- which(stations == this_station)
  if(progind %% 10 == 0){
    cat(progind, "of", length(stations), "started at", paste(Sys.time()), "\n")
  }
  
  ind_this_station <- which(data_train$station == this_station)
  
  # omit stations with too few observation to successfully estimate model
  if(length(ind_this_station) <= 10){next}
  
  # estimate model on training data
  data_train_thisstation <- subset(data_train, station == this_station)[!(names(data_train) %in% c("obs", "date", "station"))]
  y_train <- subset(data_train, station == this_station)$obs
    
  qrf_model <- quantile_forest(data_train_thisstation, 
                               y_train,
                               quantiles = qt_levels)
  
  # compute quantiles on evaluation data
  data_test_thisstation <- subset(data_eval, station == this_station)[!(names(data_eval) %in% c("obs", "date", "station"))]
    
  grF_prediction <-   predict(grF_model,
                              data_test_thisstation,
                              quantiles = qt_levels)
  
  ind_this_station_eval <- which(data_eval$station == this_station)
  
  qts_save[ind_this_station_eval,] <- grF_prediction
}


ind_noNAinRow <- which(!apply(qts_save, 1, anyNA))
ind_use <- intersect(which(!is.na(data_eval$obs)), ind_noNAinRow)

dim(qts_save)[1] - length(ind_use)

grf_crps <- crps_sample(y = data_eval$obs[ind_use],
                        dat = qts_save[ind_use,])

summary(grf_crps) # mean 1.2233  
