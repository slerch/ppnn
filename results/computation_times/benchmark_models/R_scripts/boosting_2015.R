cat(paste("Now running: Boosting (training 2015), started", Sys.time(), "\n"))
cat(paste("... should take around 15 min", "\n"))


train_end <- as.Date("2016-01-01 00:00 UTC") - days(2)
train_start <- as.Date("2015-01-01 00:00 UTC") 
data_train_all <- subset(data, date >= train_start & date <= train_end)

eval_start <- as.Date("2016-01-01 00:00 UTC")
eval_end <- as.Date("2016-12-31 00:00 UTC")
eval_dates <- seq(eval_start, eval_end, by = "1 day")

data_eval_all <- subset(data, date >= eval_start & date <= eval_end)

out_loc <- rep(NA, nrow(data_eval_all))
out_sc <- rep(NA, nrow(data_eval_all))

stations_list <- unique(data$station)

t1 <- Sys.time()

for(this_station in stations_list){
  # progress indicator
  # progind <- which(stations_list == this_station)
  # if(progind %% 10 == 0){
  #   cat(progind, "of", length(stations_list), "started at", paste(Sys.time()), "\n")
  # }
  
  # data_train <- subset(data, date >= train_start & date <= train_end & station == this_station)
  data_train <- subset(data_train_all, station == this_station)
  
  # remove incomplete cases (= NA obs or fc)
  data_train <- data_train[complete.cases(data_train), ]
  if(nrow(data_train) < 10){
    next
  }
  
  ## NOTE: boosting is only implemented for link.scale = "log", otherwise cryptic error message 
  crch_model <- crch(obs ~ .|.,
                     data = data_train[,-which(names(data) %in% c("date", "station"))],
                     dist = "gaussian",
                     link.scale = "log",
                     method = "boosting",
                     maxit = 1000,
                     nu = 0.05,
                     mstop = "aic")
  
  data_eval <- subset(data_eval_all, station == this_station)
  
  loc <- as.numeric(predict(crch_model, data_eval, type = "location"))
  sc <- as.numeric(predict(crch_model, data_eval, type = "scale"))
  
  # save parameters
  ind_this_station <- which(data_eval_all$station == this_station)
  
  out_loc[ind_this_station] <- loc 
  out_sc[ind_this_station] <- sc
}

t2 <- Sys.time()
time_elapsed <- difftime(t2, t1, units='mins')
cat("... finished, took", round(time_elapsed,2), "minutes", "\n")

