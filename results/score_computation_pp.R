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

NA_cases_training15 <- which(is.na(df_crps$emos_loc_15))

## remove remaining NA cases in qrf and local model
df_crps_withNAs <- df_crps
df_crps <- df_crps_withNAs[complete.cases(df_crps_withNAs),]

apply(df_crps[,4:ncol(df_crps)], 2, function(x) sum(is.na(x)))
apply(df_crps[,4:ncol(df_crps)], 2, mean)
round(apply(df_crps[,4:ncol(df_crps)], 2, mean),2)

df_crps_15 <- df_crps

## ---------- AE computation ---------- ##

ae_norm <- function(y, m, s){
  abs(y-m)
}
ae_sample <- function(y, dat){
  mm <- apply(dat, 1, median)
  return(abs(y-mm))
}

df_ae <- data_eval

# network models
df_ae$fc_15 <- ae_norm(data_eval$obs, fc_15$mean, fc_15$std)
df_ae$fc_aux_15 <- ae_norm(data_eval$obs, fc_aux_15$mean, fc_aux_15$std)
df_ae$fc_emb_15 <- ae_norm(data_eval$obs, fc_emb_15$mean, fc_emb_15$std)
df_ae$fc_aux_emb_15 <- ae_norm(data_eval$obs, fc_aux_emb_15$mean, fc_aux_emb_15$std)
df_ae$nn_aux_emb_15 <- ae_norm(data_eval$obs, nn_aux_emb_15$mean, nn_aux_emb_15$std)

apply(df_ae[,4:ncol(df_ae)], 2, mean, na.rm = TRUE)
apply(df_ae[,4:ncol(df_ae)], 2, function(x) sum(is.na(x)))

# benchmark models
df_ae$emos_gl_15 <- ae_norm(data_eval$obs, emos_gl_15$mean, emos_gl_15$std)
df_ae$emos_loc_15 <- ae_norm(data_eval$obs, emos_loc_15$mean, emos_loc_15$std)
df_ae$bst_15 <- ae_norm(data_eval$obs, bst_15$mean, bst_15$std)
# qrf model: special role due to quantiles as predictions
df_ae$qrf_15 <- NA
qrf_15_matrix <- as.matrix(qrf_15[,3:ncol(qrf_15)])
use_row_qrf <- !apply(qrf_15_matrix, 1, anyNA)
tmp <- ae_sample(data_eval$obs[use_row_qrf], qrf_15_matrix[use_row_qrf,])
df_ae$qrf_15[use_row_qrf] <- tmp

apply(df_ae[,4:ncol(df_ae)], 2, mean, na.rm = TRUE)
apply(df_ae[,4:ncol(df_ae)], 2, function(x) sum(is.na(x)))

## remove remaining NA cases in qrf and local model
df_ae_withNAs <- df_ae
df_ae <- df_ae_withNAs[complete.cases(df_ae_withNAs),]

apply(df_ae[,4:ncol(df_ae)], 2, function(x) sum(is.na(x)))
apply(df_ae[,4:ncol(df_ae)], 2, mean)
round(apply(df_ae[,4:ncol(df_ae)], 2, mean),2)

df_ae_15 <- df_ae

## ---------- PIT computation ---------- ##

df_pit <- data_eval

# network models
df_pit$fc_15 <- pnorm(data_eval$obs, fc_15$mean, fc_15$std)
df_pit$fc_aux_15 <- pnorm(data_eval$obs, fc_aux_15$mean, fc_aux_15$std)
df_pit$fc_emb_15 <- pnorm(data_eval$obs, fc_emb_15$mean, fc_emb_15$std)
df_pit$fc_aux_emb_15 <- pnorm(data_eval$obs, fc_aux_emb_15$mean, fc_aux_emb_15$std)
df_pit$nn_aux_emb_15 <- pnorm(data_eval$obs, nn_aux_emb_15$mean, nn_aux_emb_15$std)

# benchmark models
df_pit$emos_gl_15 <- pnorm(data_eval$obs, emos_gl_15$mean, emos_gl_15$std)
df_pit$emos_loc_15 <- pnorm(data_eval$obs, emos_loc_15$mean, emos_loc_15$std)
df_pit$bst_15 <- pnorm(data_eval$obs, bst_15$mean, bst_15$std)
# qrf model: special role due to quantiles as predictions
df_pit$qrf_15 <- NA
qrf_15_matrix <- as.matrix(qrf_15[,3:ncol(qrf_15)])
use_row_qrf <- !apply(qrf_15_matrix, 1, anyNA)
y_use <- data_eval$obs[use_row_qrf]
pred_use <- qrf_15_matrix[use_row_qrf,]
tmp <- sapply(seq_along(1:sum(use_row_qrf)),
              function(i) rank(c(y_use[i], pred_use[i,]))[1])
df_pit$qrf_15[use_row_qrf] <- tmp

apply(df_pit[,4:ncol(df_pit)], 2, function(x) sum(is.na(x)))

## remove remaining NA cases in qrf and local model
df_pit_withNAs <- df_pit
df_pit <- df_pit_withNAs[complete.cases(df_pit_withNAs),]

apply(df_pit[,4:ncol(df_pit)], 2, function(x) sum(is.na(x)))

df_pit_15 <- df_pit

# par(mfrow=c(3,3))
# for(col in 4:12){
#   if(col < 12){
#     hist(df_pit[,col], freq = FALSE, main = names(df_pit[col]), ylim = c(0,1.8))
#     abline(h = 1, lty = 2)
#   }
#   if(col == 12){
#     hist(df_pit[,col], freq = FALSE, main = names(df_pit[col]))
#     abline(h = 1/51, lty = 2)
#   }
# }

## ------------------------------------- ##
## 2007-2015 training
## ------------------------------------- ##

rm(list=setdiff(ls(), c("dir", "data_eval", "data_eval_all", "NA_cases_training15",
                        "df_crps_15", "df_pit_15", "ind_use", "csv_dir", "df_ae_15",
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

df_crps <- data_eval

# network models
df_crps$fc_0715 <- crps_norm(data_eval$obs, fc_0715$mean, fc_0715$std)
df_crps$fc_aux_0715 <- crps_norm(data_eval$obs, fc_aux_0715$mean, fc_aux_0715$std)
df_crps$fc_emb_0715 <- crps_norm(data_eval$obs, fc_emb_0715$mean, fc_emb_0715$std)
df_crps$fc_aux_emb_0715 <- crps_norm(data_eval$obs, fc_aux_emb_0715$mean, fc_aux_emb_0715$std)
df_crps$nn_aux <- crps_norm(data_eval$obs, nn_aux_0715$mean, nn_aux_0715$std)
df_crps$nn_aux_emb_0715 <- crps_norm(data_eval$obs, nn_aux_emb_0715$mean, nn_aux_emb_0715$std)

apply(df_crps[,4:ncol(df_crps)], 2, mean, na.rm = TRUE)
apply(df_crps[,4:ncol(df_crps)], 2, function(x) sum(is.na(x)))

# benchmark models
df_crps$emos_gl_0715 <- crps_norm(data_eval$obs, emos_gl_0715$mean, emos_gl_0715$std)
df_crps$emos_loc_0715 <- crps_norm(data_eval$obs, emos_loc_0715$mean, emos_loc_0715$std)
df_crps$bst_0715 <- crps_norm(data_eval$obs, bst_0715$mean, bst_0715$std)
# qrf model: special role due to quantiles as predictions
df_crps$qrf_0715 <- NA
qrf_0715_matrix <- as.matrix(qrf_0715[,3:ncol(qrf_0715)])
use_row_qrf <- !apply(qrf_0715_matrix, 1, anyNA)
tmp <- crps_sample(data_eval$obs[use_row_qrf], qrf_0715_matrix[use_row_qrf,])
df_crps$qrf_0715[use_row_qrf] <- tmp

apply(df_crps[,4:ncol(df_crps)], 2, mean, na.rm = TRUE)
apply(df_crps[,4:ncol(df_crps)], 2, function(x) sum(is.na(x)))

## only 367 missing cases here, but 551 before; remove those 551 from this data set as well
## to ensure that evaluation takes place on the same data
NA_cases_training0715 <- which(is.na(df_crps$emos_loc_0715))
any(!(NA_cases_training0715 %in% NA_cases_training15))

noNA_cases_training15 <- which(!(1:nrow(df_crps) %in% NA_cases_training15)) 

## remove remaining NA cases in qrf and local model
df_crps_withNAs <- df_crps
df_crps <- df_crps_withNAs[noNA_cases_training15,]

apply(df_crps[,4:ncol(df_crps)], 2, function(x) sum(is.na(x)))
apply(df_crps[,4:ncol(df_crps)], 2, mean)
round(apply(df_crps[,4:ncol(df_crps)], 2, mean),2)

df_crps_0715 <- df_crps


## ---------- AE computation ---------- ##

ae_norm <- function(y, m, s){
  abs(y-m)
}
ae_sample <- function(y, dat){
  mm <- apply(dat, 1, median)
  return(abs(y-mm))
}

df_ae <- data_eval

# network models
df_ae$fc_0715 <- ae_norm(data_eval$obs, fc_0715$mean, fc_0715$std)
df_ae$fc_aux_0715 <- ae_norm(data_eval$obs, fc_aux_0715$mean, fc_aux_0715$std)
df_ae$fc_emb_0715 <- ae_norm(data_eval$obs, fc_emb_0715$mean, fc_emb_0715$std)
df_ae$fc_aux_emb_0715 <- ae_norm(data_eval$obs, fc_aux_emb_0715$mean, fc_aux_emb_0715$std)
df_ae$nn_aux_0715 <- ae_norm(data_eval$obs, nn_aux_0715$mean, nn_aux_0715$std)
df_ae$nn_aux_emb_0715 <- ae_norm(data_eval$obs, nn_aux_emb_0715$mean, nn_aux_emb_0715$std)

apply(df_ae[,4:ncol(df_ae)], 2, mean, na.rm = TRUE)
apply(df_ae[,4:ncol(df_ae)], 2, function(x) sum(is.na(x)))

# benchmark models
df_ae$emos_gl_0715 <- ae_norm(data_eval$obs, emos_gl_0715$mean, emos_gl_0715$std)
df_ae$emos_loc_0715 <- ae_norm(data_eval$obs, emos_loc_0715$mean, emos_loc_0715$std)
df_ae$bst_0715 <- ae_norm(data_eval$obs, bst_0715$mean, bst_0715$std)
# qrf model: special role due to quantiles as predictions
df_ae$qrf_0715 <- NA
qrf_0715_matrix <- as.matrix(qrf_0715[,3:ncol(qrf_0715)])
use_row_qrf <- !apply(qrf_0715_matrix, 1, anyNA)
tmp <- ae_sample(data_eval$obs[use_row_qrf], qrf_0715_matrix[use_row_qrf,])
df_ae$qrf_0715[use_row_qrf] <- tmp

apply(df_ae[,4:ncol(df_ae)], 2, mean, na.rm = TRUE)
apply(df_ae[,4:ncol(df_ae)], 2, function(x) sum(is.na(x)))

## remove remaining NA cases in qrf and local model
df_ae_withNAs <- df_ae
df_ae <- df_ae_withNAs[noNA_cases_training15,]

apply(df_ae[,4:ncol(df_ae)], 2, function(x) sum(is.na(x)))
apply(df_ae[,4:ncol(df_ae)], 2, mean)
round(apply(df_ae[,4:ncol(df_ae)], 2, mean),2)

df_ae_0715 <- df_ae

## ---------- PIT computation ---------- ##

df_pit <- data_eval

# network models
df_pit$fc_0715 <- pnorm(data_eval$obs, fc_0715$mean, fc_0715$std)
df_pit$fc_aux_0715 <- pnorm(data_eval$obs, fc_aux_0715$mean, fc_aux_0715$std)
df_pit$fc_emb_0715 <- pnorm(data_eval$obs, fc_emb_0715$mean, fc_emb_0715$std)
df_pit$fc_aux_emb_0715 <- pnorm(data_eval$obs, fc_aux_emb_0715$mean, fc_aux_emb_0715$std)
df_pit$nn_aux_0715 <- pnorm(data_eval$obs, nn_aux_0715$mean, nn_aux_0715$std)
df_pit$nn_aux_emb_0715 <- pnorm(data_eval$obs, nn_aux_emb_0715$mean, nn_aux_emb_0715$std)

# benchmark models
df_pit$emos_gl_0715 <- pnorm(data_eval$obs, emos_gl_0715$mean, emos_gl_0715$std)
df_pit$emos_loc_0715 <- pnorm(data_eval$obs, emos_loc_0715$mean, emos_loc_0715$std)
df_pit$bst_0715 <- pnorm(data_eval$obs, bst_0715$mean, bst_0715$std)
# qrf model: special role due to quantiles as predictions
df_pit$qrf_0715 <- NA
qrf_0715_matrix <- as.matrix(qrf_0715[,3:ncol(qrf_0715)])
use_row_qrf <- !apply(qrf_0715_matrix, 1, anyNA)
y_use <- data_eval$obs[use_row_qrf]
pred_use <- qrf_0715_matrix[use_row_qrf,]
tmp <- sapply(seq_along(1:sum(use_row_qrf)),
              function(i) rank(c(y_use[i], pred_use[i,]))[1])
df_pit$qrf_0715[use_row_qrf] <- tmp

apply(df_pit[,4:ncol(df_pit)], 2, function(x) sum(is.na(x)))

## remove remaining NA cases in qrf and local model
df_pit_withNAs <- df_pit
df_pit <- df_pit_withNAs[noNA_cases_training15,]

apply(df_pit[,4:ncol(df_pit)], 2, function(x) sum(is.na(x)))

df_pit_0715 <- df_pit

# par(mfrow=c(3,3))
# for(col in 4:12){
#   if(col < 12){
#     hist(df_pit[,col], freq = FALSE, main = names(df_pit[col]), ylim = c(0,1.8))
#     abline(h = 1, lty = 2)
#   }
#   if(col == 12){
#     hist(df_pit[,col], freq = FALSE, main = names(df_pit[col]))
#     abline(h = 1/51, lty = 2)
#   }
# }

save(df_crps_15, df_crps_0715, df_pit_15, df_pit_0715, df_ae_15, df_ae_0715,
     ind_use, noNA_cases_training15,
     file = "CRPS_AE_PIT_all.Rdata")
