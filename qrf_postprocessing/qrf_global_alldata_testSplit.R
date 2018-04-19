## Implementation of quantile regression forest approach to post-processing
## proposed by Taillardat et al (2016, Monthly Weather Review)
## partly based on R code provided by Maxime Taillardat

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

# ## splitting into a specific number of parts here
# ## however, this is less straightforward compared to test how large indiv parts can be, see below
# ## split data_train into 10 seperate RF models, and combine them later on
# ## this could also be done randomly by permuting the indices first
# permuted_rowIndices <- sample(1:nrow(data_train), size = nrow(data_train), replace = FALSE)
# nparts <- 100
# # following row would split into parts of size npart, which is not what we want
# # rowInd_splitted <- split(permuted_rowIndices, ceiling(seq_along(permuted_rowIndices)/nparts))
# # instead:
# part_size <- floor(length(permuted_rowIndices)/nparts)
# rowInd_splitted <- split(permuted_rowIndices, ceiling(seq_along(permuted_rowIndices)/part_size))
# str(rowInd_splitted)
# # drop last list entry so that other ones have the same size 
# rowInd_splitted[[nparts+1]] <- NULL
# str(rowInd_splitted)

## split data_train into parts of size say 5000
## this could also be done randomly by permuting the indices first
permuted_rowIndices <- sample(1:nrow(data_train), size = nrow(data_train), replace = FALSE)
part_size <- 5000 
rowInd_splitted <- split(permuted_rowIndices, ceiling(seq_along(permuted_rowIndices)/part_size))
str(rowInd_splitted)
# drop last list entry so that other ones have the same size 
rowInd_splitted[[length(rowInd_splitted)]] <- NULL
str(rowInd_splitted)


# now, we have a list of indices to train separate qrF models on
rF_list <- list()
nparts <- length(rowInd_splitted)
for(ID in 1:nparts){
  
  cat(ID, "of", nparts, "started at", paste(Sys.time()), "\n")
  rows_use <- rowInd_splitted[[ID]]
  data_train_use <- data_train[rows_use,]
  
  qrF_model <- quantregForest(x = data_train_use[!(names(data_train) %in% c("date", "station", "obs"))], 
                              y = data_train_use$obs,
                              ntree = 1000,
                              nodesize = 5,
                              nthreads = 1,
                              mtry = 20,
                              replace = TRUE)
  
  rF_list[[ID]] <- qrF_model
}


rF_model_comb <- do.call("combine", rF_list)
