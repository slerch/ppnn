cat(paste("Now running: QRF (training 2015), started", Sys.time(), "\n"))
cat(paste("... should take around 8 min", "\n"))


start_eval <- as.Date("2016-01-01 00:00", tz = "UTC")
end_eval <- as.Date("2016-12-31 00:00", tz = "UTC") 

data_eval_all <- subset(data, date >= start_eval & date <= end_eval)

end_train <- start_eval - hours(48)
start_train <- as.Date("2015-01-01 00:00", tz  = "UTC")

data_train_all <- subset(data, date >= start_train & date <= end_train )

# quantile levels
qt_levels <- seq(1/51, 50/51, by = 1/51)
qts_save <- matrix(NA, 
                   nrow = nrow(data_eval_all), 
                   ncol = length(qt_levels))

stations <- unique(data$station)

t1 <- Sys.time()

for(this_station in stations){

  data_train <- subset(data_train_all, station == this_station)
  
  data_train <- data_train[complete.cases(data_train),]
  
  data_eval <- subset(data_eval_all, station == this_station)
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
  

  qts_save[which(data_eval_all$station == this_station),] <- qrF_prediction
}

t2 <- Sys.time()
time_elapsed <- difftime(t2, t1, units='mins')
cat("... finished, took", round(time_elapsed,2), "minutes", "\n")
