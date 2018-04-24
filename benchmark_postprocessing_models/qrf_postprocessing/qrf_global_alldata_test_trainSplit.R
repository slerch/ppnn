## Implementation of quantile regression forest approach to post-processing
## proposed by Taillardat et al (2016, Monthly Weather Review)
## partly based on R code provided by Maxime Taillardat

##
## running global QRF model for all data requires splitting training data
## and fitting multiple models on random subsets
## However:
## bug in "combine()" function for quantregForest objects
## I have contacted the author of the package who will look for a solution
##

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

data_eval <- subset(data, date >= start_eval & date <= end_eval)

end_train <- start_eval - hours(48)
start_train <- data$date[1]

qt_levels <- seq(1/21, 20/21, by = 1/21)
qts_save <- matrix(NA, 
                   nrow = nrow(data_eval), 
                   ncol = length(qt_levels))

data_train <- subset(data, date >= start_train & date <= end_train)
data_train <- data_train[complete.cases(data_train),]

# split up data set into smaller chunks
# and fit rF model on each chunk individually
# this avoids lack of memory issues when fitting models with large data sets
# also allows straightforward parallelization of experiments on cluster

## split data_train into parts of size say 5000
## this could also be done randomly by permuting the indices first
permuted_rowIndices <- sample(1:nrow(data_train), size = nrow(data_train), replace = FALSE)
part_size <- 5000 
rowInd_splitted <- split(permuted_rowIndices, ceiling(seq_along(permuted_rowIndices)/part_size))
str(rowInd_splitted)
# drop last list entry so that other ones have the same size 
rowInd_splitted[[length(rowInd_splitted)]] <- NULL
str(rowInd_splitted)


# # now, we have a list of indices to train separate qrF models on
# rF_list <- list()
# nparts <- length(rowInd_splitted)
# for(ID in 1:nparts){
#   
#   cat(ID, "of", nparts, "started at", paste(Sys.time()), "\n")
#   rows_use <- rowInd_splitted[[ID]]
#   data_train_use <- data_train[rows_use,]
#   
#   qrF_model <- quantregForest(x = data_train_use[!(names(data_train) %in% c("date", "station", "obs"))], 
#                               y = data_train_use$obs,
#                               ntree = 1000,
#                               nodesize = 5,
#                               nthreads = 1,
#                               mtry = 20,
#                               replace = TRUE)
#   
#   rF_list[[ID]] <- qrF_model
# }
# 
# 
# rF_model_comb <- do.call("combine", rF_list)

# better:

nparts <- length(rowInd_splitted)

fit_model <- function(ID){
  cat(ID, "of", nparts, "started at", paste(Sys.time()), "\n")
  rows_use <- rowInd_splitted[[ID]]
  data_train_use <- data_train[rows_use,]
  
  qrF_model <- quantregForest(x = data_train_use[!(names(data_train) %in% c("date", "station", "obs"))], 
                              y = data_train_use$obs,
                              ntree = 250,
                              nodesize = 5,
                              nthreads = 1,
                              mtry = 20,
                              replace = TRUE)
  
  return(qrF_model)
}

rF_list <- lapply(as.list(1:nparts), fit_model)

rF_model_comb <- do.call("combine", rF_list)

## should work like this, however, combined model results in strangely sharp predictions,
## or even only 0's for all quantiles

data_eval <- subset(data, date >= start_eval & date <= end_eval)
data_eval <- data_eval[1:500,]

nn <- 20
qt_levels <- seq(1/(nn+1), nn/(nn+1), by = 1/(nn+1)) # quantile levels

qrF_prediction <-   predict(rF_model_comb,
                            newdata = data_eval[, !(names(data_eval) %in% c("date", "station", "obs"))],
                            what = qt_levels,
                            all = TRUE)

head(qrF_prediction)
