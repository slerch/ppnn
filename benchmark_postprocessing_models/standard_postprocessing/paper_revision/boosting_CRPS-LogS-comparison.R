## Data preparation: Run

# rm(list=ls())
# 
# data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
# load(paste0(data_dir, "data_all.Rdata"))
# 
# # remove sm data (missing values)
# data$sm_mean <- NULL
# data$sm_var <- NULL
# head(data)
# 
# library(scoringRules)
# library(lubridate)
# library(crch)
# 
# train_end <- as.Date("2016-01-01 00:00 UTC") - days(2)
# train_start <- data$date[1]
# 
# data_train_all <- subset(data, date >= train_start & date <= train_end)
# 
# eval_start <- as.Date("2016-01-01 00:00 UTC")
# eval_end <- as.Date("2016-12-31 00:00 UTC")
# eval_dates <- seq(eval_start, eval_end, by = "1 day")
# 
# data_eval_all <- subset(data, date >= eval_start & date <= eval_end)
# 
# out_loc <- rep(NA, nrow(data_eval_all))
# out_sc <- rep(NA, nrow(data_eval_all))
# 
# stations_list <- unique(data$station)

## Logs: run 

# ncoef_loc <- NULL
# ncoef_sc <- NULL
# 
# for(this_station in stations_list){
#   # progress indicator
#   progind <- which(stations_list == this_station)
#   if(progind %% 10 == 0){
#     cat(progind, "of", length(stations_list), "started at", paste(Sys.time()), "\n")
#   }
#   
#   # data_train <- subset(data, date >= train_start & date <= train_end & station == this_station)
#   data_train <- subset(data_train_all, station == this_station)
#   
#   # remove incomplete cases (= NA obs or fc)
#   data_train <- data_train[complete.cases(data_train), ]
#   if(nrow(data_train) < 10){
#     next
#   }
#   
#   ## NOTE: boosting is only implemented for link.scale = "log", otherwise cryptic error message 
#   crch_model <- crch(obs ~ .|.,
#                      data = data_train[,-which(names(data) %in% c("date", "station"))],
#                      dist = "gaussian",
#                      link.scale = "log",
#                      method = "boosting",
#                      maxit = 1000,
#                      nu = 0.05,
#                      mstop = "aic")
#   
#   ncoef_loc[progind] <- sum(crch_model$coefficients$location > 0)
#   ncoef_sc[progind] <- sum(crch_model$coefficients$scale > 0)
# }
# 
# save(ncoef_loc, ncoef_sc, file = "ncoef_LogS.Rdata")

## CRPS: run

# ncoef_loc <- NULL
# ncoef_sc <- NULL
# 
# for(this_station in stations_list){
#   # progress indicator
#   progind <- which(stations_list == this_station)
#   if(progind %% 10 == 0){
#     cat(progind, "of", length(stations_list), "started at", paste(Sys.time()), "\n")
#   }
#   
#   # data_train <- subset(data, date >= train_start & date <= train_end & station == this_station)
#   data_train <- subset(data_train_all, station == this_station)
#   
#   # remove incomplete cases (= NA obs or fc)
#   data_train <- data_train[complete.cases(data_train), ]
#   if(nrow(data_train) < 10){
#     next
#   }
#   
#   ## NOTE: boosting is only implemented for link.scale = "log", otherwise cryptic error message 
#   crch_model <- crch(obs ~ .|.,
#                      data = data_train[,-which(names(data) %in% c("date", "station"))],
#                      dist = "gaussian",
#                      link.scale = "log",
#                      method = "boosting",
#                      type = "crps",
#                      maxit = 1000,
#                      nu = 0.05,
#                      mstop = "aic")
#   
#   ncoef_loc[progind] <- sum(crch_model$coefficients$location > 0)
#   ncoef_sc[progind] <- sum(crch_model$coefficients$scale > 0)
# }
# save(ncoef_loc, ncoef_sc, file = "ncoef_CRPS.Rdata")

## Plot
rm(list=ls())

load("ncoef_CRPS.Rdata")
ncoef_loc_crps <- ncoef_loc
ncoef_sc_crps <- ncoef_sc

load("ncoef_LogS.Rdata")
ncoef_loc_logs <- ncoef_loc
ncoef_sc_logs <- ncoef_sc

pdf("ncoef_boosting_CRPS-LogS.pdf", width = 10, height = 5, pointsize = 12)

par(mfrow=c(1,2))

# plot for location parameter
ploc_logs <- hist(ncoef_loc_logs, breaks = seq(0,25,1), plot = FALSE)                    
ploc_crps <- hist(ncoef_loc_crps, breaks = seq(0,25,1), plot = FALSE)                     
plot(ploc_logs, col=rgb(0,0,1,1/4), xlim=c(0,25), ylim = c(0,150), 
      main = "Location", xlab = "Number of selected predictors")  # first histogram
plot(ploc_crps, col=rgb(1,0,0,1/4), xlim=c(0,25), add=T)  # second
legend("topright", c("LogS", "CRPS"), fill=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)))

# plot for scale parameter
psc_logs <- hist(ncoef_sc_logs, breaks = seq(0,25,1), plot = FALSE)                    
psc_crps <- hist(ncoef_sc_crps, breaks = seq(0,25,1), plot = FALSE)                     
plot(psc_logs, col=rgb(0,0,1,1/4), xlim=c(0,25), ylim = c(0,150), 
     main = "Scale", xlab = "Number of selected predictors")  # first histogram
plot(psc_crps, col=rgb(1,0,0,1/4), xlim=c(0,25), add=T)  # second
legend("topright", c("LogS", "CRPS"), fill=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)))

dev.off()