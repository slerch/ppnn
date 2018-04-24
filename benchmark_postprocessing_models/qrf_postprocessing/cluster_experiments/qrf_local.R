## implementation on cluster

# re-write as function of tuning parameters,
# and run over expand.grid of those values with as SGE tasks
# save output to Rdata files (mean crps will suffice if files are large)

rm(list=ls())

data_dir <- "..."
load(paste0(data_dir, "data_all.Rdata"))

## ignore sm_mean and sm_var due to missing and NA values leading to problems in estimation
data$sm_mean <- NULL
data$sm_var <- NULL

library(randomForest, lib = "...")
library(RColorBrewer, lib = "...")
library(quantregForest, lib = "...")
library(scoringRules, lib = "...")
library(lubridate, lib = "...")

# train and eval dates (any date in date = valid date of fc and!)
start_eval <- as.Date("2016-01-01 00:00", tz = "UTC")
end_eval <- as.Date("2016-12-31 00:00", tz = "UTC") 

data_eval_all <- subset(data, date >= start_eval & date <= end_eval)

end_train <- start_eval - hours(48)
start_train <- as.Date("2015-01-01 00:00", tz = "UTC")

# quantile levels
qt_levels <- seq(1/51, 50/51, by = 1/51)
qts_save <- matrix(NA, 
                   nrow = nrow(data_eval_all), 
                   ncol = length(qt_levels))

stations <- unique(data$station)

# see below for description

fit_qrf <- function(ntree, nodesize, mtry, replace){
  
  Routname1 <- paste0("qrf_local_ntree", ntree, "_nodesize", nodesize, "_mtry", mtry, "_repl", replace, ".Rout")
  Routname <- paste0("...", Routname1)
  
  sink(Routname)
  
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
                                nodesize = nodesize,
                                ntree = ntree,
                                mtry = mtry,
                                replace = replace)
    
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
  
  savename <- paste0("qrf_local_ntree", ntree, "_nodesize", nodesize, "_mtry", mtry, "_repl", replace, ".Rdata")
  fname <- paste0("...", savename)
  
  save(qrf_crps, file = fname)
  
  sink()

}

# ntree: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.
#     default: 500
# nodesize: Minimum size of terminal nodes
#     default: 5
# mtry: Number of variables randomly sampled as candidates at each split
#     default: max(floor(ncol(x)/3), 1) = 11
# replace: Should sampling of cases be done with or without replacement?
#     default: 1 (or 0)

ntree.try <- c(125,250,500,1000,2500)
nodesize.try <- c(5,10,15,20,50)
mtry.try <- c(10,15,20,25)
replace.try <- c(0,1)

pars <- expand.grid(ntree.try, nodesize.try, mtry.try, replace.try)

nrow(pars) # --> 200 tasks

ID <- as.numeric(as.character(Sys.getenv(c("SGE_TASK_ID"))))

fit_qrf(ntree = pars[ID,1], nodesize = pars[ID,2], mtry = pars[ID,3], replace = pars[ID,4])
