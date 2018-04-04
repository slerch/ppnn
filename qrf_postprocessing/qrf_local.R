## Implementation of quantile regression forest approach to post-processing
## proposed by Taillardat et al (2016, Monthly Weather Review)
## partly based on R code provided by Maxime Taillardat
## here: local = one model for every station

rm(list=ls())

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
  data_train_thisstation <- subset(data_train, station == this_station)
  
  qrF_model <- quantregForest(x = data_train_thisstation[,2:(ncol(data_train)-2)], 
                              y = data_train_thisstation$obs,
                              mtry = ceiling(ncol(data_train_thisstation[,2:(ncol(data_train)-2)])/3),
                              nodesize = 5,
                              ntree = 250,
                              importance = FALSE)
  
  # compute quantiles on evaluation data
  
  ind_this_station_eval <- which(data_eval$station == this_station)
  
  qrF_prediction <-   predict(qrF_model,
                              newdata = subset(data_eval, station == this_station),
                              what = qt_levels,
                              all = TRUE)
  
  qts_save[ind_this_station_eval,] <- qrF_prediction
}


ind_noNAinRow <- which(!apply(qts_save, 1, anyNA))
ind_use <- intersect(which(!is.na(data_eval$obs)), ind_noNAinRow)

dim(qts_save)[1] - length(ind_use)

qrf_crps <- crps_sample(y = data_eval$obs[ind_use],
                        dat = qts_save[ind_use,])

summary(qrf_crps)

# nodesize = 10, ntree = 250, --> mean CRPS 0.9964
# nodesize = 5, ntree = 250, --> mean CRPS 0.9906
# nodesize = 5, maxnodes = 20, ntree = 250, --> mean CRPS 1.058
# all to default: --> mean CRPS 0.9972