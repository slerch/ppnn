rm(list=ls())

library(scoringRules)

csv_dir <- "/home/sebastian/Projects/PP_NN/code/results/csv_files/"
names_keep <- c("date", "station_id", "mean", "std")

# load network model results

fc_15 <- read.csv(paste0(csv_dir, "fc_15.csv"))
# drop un-used names
fc_15[,!(names(fc_15) %in% names_keep)] <- NULL

fc_aux_15 <- read.csv(paste0(csv_dir, "fc_aux_15.csv"))
# drop un-used names
fc_aux_15[,!(names(fc_aux_15) %in% names_keep)] <- NULL

fc_emb_15 <- read.csv(paste0(csv_dir, "fc_emb_15.csv"))
# drop un-used names
fc_emb_15[,!(names(fc_emb_15) %in% names_keep)] <- NULL

fc_aux_emb_15 <- read.csv(paste0(csv_dir, "fc_aux_emb_15.csv"))
# drop un-used names
fc_aux_emb_15[,!(names(fc_aux_emb_15) %in% names_keep)] <- NULL

nn_aux_emb_15 <- read.csv(paste0(csv_dir, "nn_aux_emb_15.csv"))
# drop un-used names
nn_aux_emb_15[,!(names(nn_aux_emb_15) %in% names_keep)] <- NULL

## benchmark model data

library(lubridate)

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_all.Rdata"))
start_eval <- as.Date("2016-01-01 00:00", tz = "UTC")
end_eval <- as.Date("2016-12-31 00:00", tz = "UTC") 
data_eval_all <- subset(data, date >= start_eval & date <= end_eval)
data_eval_all[, !(names(data_eval_all) %in% c("date", "station", "obs"))] <- NULL
rm(data)

ind_use <- which(!is.na(data_eval$obs))
data_eval <- data_eval_all[ind_use,]


