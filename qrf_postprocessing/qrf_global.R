## Implementation of quantile regression forest approach to post-processing
## proposed by Taillardat et al (2016, Monthly Weather Review)
## partly based on R code provided by Maxime Taillardat

rm(list=ls())

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_prep_2015-16.Rdata"))

# remove sm_mean and sm_var due to missing values
data_train$sm_mean <- NULL
data_train$sm_var <- NULL
data_eval$sm_mean <- NULL
data_eval$sm_var <- NULL

library(quantregForest)
library(scoringRules)

# tuning parameters based on code by Maxime Taillardat
# estimation takes ~ 3 hours on my laptop
qrF_model <- quantregForest(x = data_train[,2:(ncol(data_train)-2)], y = data_train$obs,
                            mtry = ceiling(ncol(data_train[,2:(ncol(data_train)-2)])/3),
                            nodesize = 10,
                            ntree = 250,
                            importance = TRUE)
  
qt_levels <- seq(1/51, 50/51, by = 1/51)
qts_save <- matrix(NA, nrow = nrow(data_eval), ncol = length(qt_levels))

qrF_prediction <-   predict(qrF_model,
                            data_eval,
                            what = qt_levels,
                            all = TRUE)

ind_noNAinRow <- which(!apply(qrF_prediction, 1, anyNA))
ind_use <- intersect(which(!is.na(data_eval$obs)), ind_noNAinRow)

qrf_crps <- crps_sample(y = data_eval$obs[ind_use],
                        dat = qrF_prediction[ind_use,])

summary(qrf_crps)

## with nodesize = 10, ntree = 250; mean 0.9890

qrF_model$importance
importance(qrF_model, type=1, scale=TRUE)


save(qrF_model, file = "qrF_model.Rdata")
