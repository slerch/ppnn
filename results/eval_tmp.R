rm(list=ls())

## ---------- load data ---------- ##

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

ind_use <- which(!is.na(data_eval_all$obs))
data_eval <- data_eval_all[ind_use,]


# benchmark models

Rdata_dir <- "/home/sebastian/Projects/PP_NN/code/results/Rdata_files/"

load(paste0(Rdata_dir, "emos_global_fixed_2015.Rdata"))
emos_gl_15 <- emos_global_fixed_2015[ind_use,]
rm(emos_global_fixed_2015)

load(paste0(Rdata_dir, "emos_global_window_2015.Rdata"))
emos_gl_w_15 <- emos_global_window_2015[ind_use,]
rm(emos_global_window_2015)

load(paste0(Rdata_dir, "emos_local_fixed_2015.Rdata"))
emos_loc_15 <- emos_local_fixed_2015[ind_use,]
rm(emos_local_fixed_2015)

load(paste0(Rdata_dir, "emos_local_window_2015.Rdata"))
emos_loc_w_15 <- emos_local_window_2015[ind_use,]
rm(emos_local_window_2015)

load(paste0(Rdata_dir, "boosting_local_2015.Rdata"))
bst_15 <- boosting_local_2015[ind_use,]
rm(boosting_local_2015)

load(paste0(Rdata_dir, "qrf_local_2015.Rdata"))
qrf_15 <- qrf_pred_2015[ind_use,]
rm(qrf_pred_2015)

# check consistency with network model results
dim(data_eval)
dim(fc_15)
dim(emos_gl_15)
any(as.Date(fc_15$date) != data_eval$date)
any(emos_gl_15$date != data_eval$date)
any(as.Date(fc_15$date) != emos_gl_15$date)


## ---------- CRPS computation ---------- ##

library(scoringRules)

df_crps <- data_eval

# network models
df_crps$fc_15 <- crps_norm(data_eval$obs, fc_15$mean, fc_15$std)
df_crps$fc_aux_15 <- crps_norm(data_eval$obs, fc_aux_15$mean, fc_aux_15$std)
df_crps$fc_emb_15 <- crps_norm(data_eval$obs, fc_emb_15$mean, fc_emb_15$std)
df_crps$fc_aux_emb_15 <- crps_norm(data_eval$obs, fc_aux_emb_15$mean, fc_aux_emb_15$std)
df_crps$nn_aux_emb_15 <- crps_norm(data_eval$obs, nn_aux_emb_15$mean, nn_aux_emb_15$std)

apply(df_crps[,4:ncol(df_crps)], 2, mean, na.rm = TRUE)
apply(df_crps[,4:ncol(df_crps)], 2, function(x) sum(is.na(x)))

# benchmark models
df_crps$emos_gl_15 <- crps_norm(data_eval$obs, emos_gl_15$mean, emos_gl_15$std)
df_crps$emos_loc_15 <- crps_norm(data_eval$obs, emos_loc_15$mean, emos_loc_15$std)
df_crps$bst_15 <- crps_norm(data_eval$obs, bst_15$mean, bst_15$std)
# qrf model: special role due to quantiles as predictions
df_crps$qrf_15 <- NA
qrf_15_matrix <- as.matrix(qrf_15[,3:ncol(qrf_15)])
use_row_qrf <- !apply(qrf_15_matrix, 1, anyNA)
tmp <- crps_sample(data_eval$obs[use_row_qrf], qrf_15_matrix[use_row_qrf,])
df_crps$qrf_15[use_row_qrf] <- tmp

apply(df_crps[,4:ncol(df_crps)], 2, mean, na.rm = TRUE)
apply(df_crps[,4:ncol(df_crps)], 2, function(x) sum(is.na(x)))

## remove remaining NA cases in qrf and local model
df_crps_withNAs <- df_crps
df_crps <- df_crps_withNAs[complete.cases(df_crps_withNAs),]

apply(df_crps[,4:ncol(df_crps)], 2, function(x) sum(is.na(x)))
apply(df_crps[,4:ncol(df_crps)], 2, mean)

