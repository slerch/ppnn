## implementation on cluster

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


# see below for description

fit_qrf <- function(ntree, nodesize, mtry, replace){
  
  data_train <- subset(data, date >= start_train & date <= end_train)
  
  data_train <- data_train[complete.cases(data_train),]
  
  data_eval <- subset(data, date >= start_eval & date <= end_eval)
  # need to delete cases with NA obs in eval data to avoid error in qrf package
  # (even though not needed in predict function...)
  data_eval$obs <- NULL
  
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
  
  ind_noNAinRow <- which(!apply(qrF_prediction, 1, anyNA))
  ind_use <- intersect(which(!is.na(data_eval_all$obs)), ind_noNAinRow)
  
  qrf_crps <- crps_sample(y = data_eval_all$obs[ind_use],
                          dat = qrF_prediction[ind_use,])
  
  savename <- paste0("qrf_global_ntree", ntree, "_nodesize", nodesize, "_mtry", mtry, "_repl", replace, ".Rdata")
  fname <- paste0("...", savename)
  
  save(qrf_crps, file = fname)
}

# ntree: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.
#     default: 500
# nodesize: Minimum size of terminal nodes
#     default: 5
# mtry: Number of variables randomly sampled as candidates at each split
#     default: max(floor(ncol(x)/3), 1) = 11
# replace: Should sampling of cases be done with or without replacement?
#     default: 1 (or 0) [only use 1 due to results from local model!]

ntree.try <- c(125,250,500,1000,2500)
nodesize.try <- c(5,10,15,20,50)
mtry.try <- c(10,15,20,25)
replace.try <- c(1)

pars <- expand.grid(ntree.try, nodesize.try, mtry.try, replace.try)

nrow(pars) # --> 100 tasks

ID <- as.numeric(as.character(Sys.getenv(c("SGE_TASK_ID"))))

fit_qrf(ntree = pars[ID,1], nodesize = pars[ID,2], mtry = pars[ID,3], replace = pars[ID,4])