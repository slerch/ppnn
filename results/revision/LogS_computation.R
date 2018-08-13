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

nn_aux_15 <- read.csv(paste0(csv_dir, "nn_aux_15.csv"))
# drop un-used names
nn_aux_15[,!(names(nn_aux_15) %in% names_keep)] <- NULL

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

df_logs <- data_eval

# network models
df_logs$fc_15 <- logs_norm(data_eval$obs, fc_15$mean, fc_15$std)
df_logs$fc_aux_15 <- logs_norm(data_eval$obs, fc_aux_15$mean, fc_aux_15$std)
df_logs$fc_emb_15 <- logs_norm(data_eval$obs, fc_emb_15$mean, fc_emb_15$std)
df_logs$fc_aux_emb_15 <- logs_norm(data_eval$obs, fc_aux_emb_15$mean, fc_aux_emb_15$std)
df_logs$nn_aux_15 <- logs_norm(data_eval$obs, nn_aux_15$mean, nn_aux_15$std)
df_logs$nn_aux_emb_15 <- logs_norm(data_eval$obs, nn_aux_emb_15$mean, nn_aux_emb_15$std)

apply(df_logs[,4:ncol(df_logs)], 2, mean, na.rm = TRUE)
apply(df_logs[,4:ncol(df_logs)], 2, function(x) sum(is.na(x)))

# benchmark models
df_logs$emos_gl_15 <- logs_norm(data_eval$obs, emos_gl_15$mean, emos_gl_15$std)
df_logs$emos_loc_15 <- logs_norm(data_eval$obs, emos_loc_15$mean, emos_loc_15$std)
df_logs$bst_15 <- logs_norm(data_eval$obs, bst_15$mean, bst_15$std)
# QRF forecasts ommited, not obvious how to calculate LogS

apply(df_logs[,4:ncol(df_logs)], 2, mean, na.rm = TRUE)
apply(df_logs[,4:ncol(df_logs)], 2, function(x) sum(is.na(x)))

NA_cases_training15 <- which(is.na(df_logs$emos_loc_15))

## remove remaining NA cases in qrf and local model
df_logs_withNAs <- df_logs
df_logs <- df_logs_withNAs[complete.cases(df_logs_withNAs),]

apply(df_logs[,4:ncol(df_logs)], 2, function(x) sum(is.na(x)))
apply(df_logs[,4:ncol(df_logs)], 2, mean)
round(apply(df_logs[,4:ncol(df_logs)], 2, mean),2)
apply(df_logs[,4:ncol(df_logs)], 2, median)

df_logs_15 <- df_logs

# some very large LogS values for bst model, exclude those
df_logs_15[which(df_logs_15$bst_15 > 100),] <- NA
round(apply(df_logs_15[,4:ncol(df_logs_15)], 2, mean, na.rm = TRUE),2)

sort(round(apply(df_logs_15[,4:ncol(df_logs_15)], 2, mean, na.rm = TRUE),2))

## ------------------------------------- ##
## 2007-2015 training
## ------------------------------------- ##

rm(list=setdiff(ls(), c("dir", "data_eval", "data_eval_all", "NA_cases_training15",
                        "df_logs_15", "df_pit_15", "ind_use", "csv_dir", 
                        "names_keep", "Rdata_dir")))

# load network model results

fc_0715 <- read.csv(paste0(csv_dir, "fc_07-15.csv"))
# drop un-used names
fc_0715[,!(names(fc_0715) %in% names_keep)] <- NULL

fc_aux_0715 <- read.csv(paste0(csv_dir, "fc_aux_07-15.csv"))
# drop un-used names
fc_aux_0715[,!(names(fc_aux_0715) %in% names_keep)] <- NULL

fc_emb_0715 <- read.csv(paste0(csv_dir, "fc_emb_07-15.csv"))
# drop un-used names
fc_emb_0715[,!(names(fc_emb_0715) %in% names_keep)] <- NULL

fc_aux_emb_0715 <- read.csv(paste0(csv_dir, "fc_aux_emb_07-15.csv"))
# drop un-used names
fc_aux_emb_0715[,!(names(fc_aux_emb_0715) %in% names_keep)] <- NULL

nn_aux_0715 <- read.csv(paste0(csv_dir, "nn_aux_07-15.csv"))
# drop un-used names
nn_aux_0715[,!(names(nn_aux_0715) %in% names_keep)] <- NULL

nn_aux_emb_0715 <- read.csv(paste0(csv_dir, "nn_aux_emb_07-15.csv"))
# drop un-used names
nn_aux_emb_0715[,!(names(nn_aux_emb_0715) %in% names_keep)] <- NULL

# benchmark models

Rdata_dir <- "/home/sebastian/Projects/PP_NN/code/results/Rdata_files/"

load(paste0(Rdata_dir, "emos_global_fixed_2007-15.Rdata"))
emos_gl_0715 <- emos_global_fixed_alldata[ind_use,]
rm(emos_global_fixed_alldata)

load(paste0(Rdata_dir, "emos_global_window_2007-15.Rdata"))
emos_gl_w_0715 <- emos_global_window_alldata[ind_use,]
rm(emos_global_window_alldata)

load(paste0(Rdata_dir, "emos_local_fixed_2007-15.Rdata"))
emos_loc_0715 <- emos_local_fixed_alldata[ind_use,]
rm(emos_local_fixed_alldata)

load(paste0(Rdata_dir, "emos_local_window_2007-15.Rdata"))
emos_loc_w_0715 <- emos_local_window_alldata[ind_use,]
rm(emos_local_window_alldata)

load(paste0(Rdata_dir, "boosting_local_2007-15.Rdata"))
bst_0715 <- boosting_local_alldata[ind_use,]
rm(boosting_local_alldata)

load(paste0(Rdata_dir, "qrf_local_2007-15.Rdata"))
qrf_0715 <- qrf_pred_alldata[ind_use,]
rm(qrf_pred_alldata)

# check consistency with network model results
dim(data_eval)
dim(fc_0715)
dim(emos_gl_0715)
any(as.Date(fc_0715$date) != data_eval$date)
any(emos_gl_0715$date != data_eval$date)
any(as.Date(fc_0715$date) != emos_gl_0715$date)

## ---------- CRPS computation ---------- ##

library(scoringRules)

df_logs <- data_eval

# network models
df_logs$fc_0715 <- logs_norm(data_eval$obs, fc_0715$mean, fc_0715$std)
df_logs$fc_aux_0715 <- logs_norm(data_eval$obs, fc_aux_0715$mean, fc_aux_0715$std)
df_logs$fc_emb_0715 <- logs_norm(data_eval$obs, fc_emb_0715$mean, fc_emb_0715$std)
df_logs$fc_aux_emb_0715 <- logs_norm(data_eval$obs, fc_aux_emb_0715$mean, fc_aux_emb_0715$std)
df_logs$nn_aux_0715 <- logs_norm(data_eval$obs, nn_aux_0715$mean, nn_aux_0715$std)
df_logs$nn_aux_emb_0715 <- logs_norm(data_eval$obs, nn_aux_emb_0715$mean, nn_aux_emb_0715$std)

apply(df_logs[,4:ncol(df_logs)], 2, mean, na.rm = TRUE)
apply(df_logs[,4:ncol(df_logs)], 2, function(x) sum(is.na(x)))

# benchmark models
df_logs$emos_gl_0715 <- logs_norm(data_eval$obs, emos_gl_0715$mean, emos_gl_0715$std)
df_logs$emos_loc_0715 <- logs_norm(data_eval$obs, emos_loc_0715$mean, emos_loc_0715$std)
df_logs$bst_0715 <- logs_norm(data_eval$obs, bst_0715$mean, bst_0715$std)
# QRF model left out, not obvious how to compute LogS

apply(df_logs[,4:ncol(df_logs)], 2, mean, na.rm = TRUE)
apply(df_logs[,4:ncol(df_logs)], 2, function(x) sum(is.na(x)))

## only 367 missing cases here, but 551 before; remove those 551 from this data set as well
## to ensure that evaluation takes place on the same data
NA_cases_training0715 <- which(is.na(df_logs$emos_loc_0715))
any(!(NA_cases_training0715 %in% NA_cases_training15))

noNA_cases_training15 <- which(!(1:nrow(df_logs) %in% NA_cases_training15)) 

## remove remaining NA cases in qrf and local model
df_logs_withNAs <- df_logs
df_logs <- df_logs_withNAs[noNA_cases_training15,]

apply(df_logs[,4:ncol(df_logs)], 2, function(x) sum(is.na(x)))
apply(df_logs[,4:ncol(df_logs)], 2, mean)
round(apply(df_logs[,4:ncol(df_logs)], 2, mean),2)

df_logs_0715 <- df_logs


## overview of results
sort(round(apply(df_logs_0715[,4:ncol(df_logs_0715)], 2, mean),2))
df_logs_15 <- df_logs_15[complete.cases(df_logs_15),]
sort(round(apply(df_logs_15[,4:ncol(df_logs_15)], 2, mean),2))

save(df_logs_15, df_logs_0715, 
     ind_use, noNA_cases_training15,
     file = "/media/sebastian/Elements/Postproc_NN/LogS_all.Rdata")
