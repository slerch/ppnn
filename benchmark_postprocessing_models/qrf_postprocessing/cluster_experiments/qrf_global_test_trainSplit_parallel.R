rm(list=ls())

##
## splitting into multiple trees and combining does not work
## due to a bug in the corresponding package (strange "jumpy" quantiles)
##

data_dir <- "..."
load(paste0(data_dir, "data_all.Rdata"))

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

data_train <- subset(data, date >= start_train & date <= end_train)
data_train <- data_train[complete.cases(data_train),]

# split up data set into smaller chunks
# and fit rF model on each chunk individually
# this avoids lack of memory issues when fitting models with large data sets
# also allows straightforward parallelization of experiments on cluster

## split data_train into parts of size say 1000
## this could also be done randomly by permuting the indices first
permuted_rowIndices <- sample(1:nrow(data_train), size = nrow(data_train), replace = FALSE)
part_size <- 5000 
rowInd_splitted <- split(permuted_rowIndices, ceiling(seq_along(permuted_rowIndices)/part_size))
# drop last list entry so that other ones have the same size 
rowInd_splitted[[length(rowInd_splitted)]] <- NULL
str(rowInd_splitted)

# rF_list <- list()
nparts <- length(rowInd_splitted)
ID_list <- as.list(1:nparts)

fit_model <- function(ID){
  library(randomForest, lib = "/home/lerchsn/Rpackages/")
  library(RColorBrewer, lib = "/home/lerchsn/Rpackages/")
  library(quantregForest, lib = "/home/lerchsn/Rpackages/")
  
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

library(snow, lib.loc = "/home/lerchsn/Rpackages/")
library(Rmpi, lib.loc = "/home/lerchsn/Rpackages/")
cl <- makeCluster(12, "MPI")
clusterEvalQ(cl, Sys.info()['nodename'])
clusterExport(cl=cl, list = ls())

rF_list <- parLapply(cl = cl, x = ID_list, fun = fit_model)

rF_model_comb <- do.call("combine", rF_list)

rm(rF_list)
print(rF_model_comb)

# predict
data_eval <- subset(data, date >= start_eval & date <= end_eval)
data_eval <- data_eval[1:1000,]

qt_levels <- seq(1/21, 20/21, by = 1/21) # quantile levels

qrF_prediction <-   predict(rF_model_comb,
                            newdata = data_eval[, !(names(data_eval) %in% c("date", "station", "obs"))],
                            what = qt_levels,
                            all = TRUE)

head(qrF_prediction)


stopCluster(cl)