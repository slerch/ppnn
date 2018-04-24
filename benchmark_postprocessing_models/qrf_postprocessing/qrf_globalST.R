## Implementation of quantile regression forest approach to post-processing
## proposed by Taillardat et al (2016, Monthly Weather Review)
## partly based on R code provided by Maxime Taillardat
## globalST variant = added station information to global data

rm(list=ls())

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_all.Rdata"))
load(paste0(data_dir, "data_all_stationInfo.Rdata"))

data$station_alt <- station_alt
data$station_orog <- station_orog
data$station_lsm <- station_lsm
rm(station_alt, station_orog, station_lsm)

## ignore sm_mean and sm_var due to missing and NA values leading to problems in estimation
data$sm_mean <- NULL
data$sm_var <- NULL

library(quantregForest)
library(scoringRules)
library(lubridate)

start_eval <- as.Date("2016-01-01 00:00", tz = "UTC")
end_eval <- as.Date("2016-12-31 00:00", tz = "UTC") 

data_eval_all <- subset(data, date >= start_eval & date <= end_eval)

end_train <- start_eval - hours(48)
start_train <- as.Date("2015-01-01 00:00", tz = "UTC")

qt_levels <- seq(1/51, 50/51, by = 1/51)
qts_save <- matrix(NA, 
                   nrow = nrow(data_eval_all), 
                   ncol = length(qt_levels))

data_train <- subset(data, date >= start_train & date <= end_train)
data_train <- data_train[complete.cases(data_train),]

data_eval <- subset(data, date >= start_eval & date <= end_eval)

# model
qrF_model <- quantregForest(x = data_train[, !(names(data_train) %in% c("date", "station", "obs"))], 
                            y = data_train$obs)

# compute quantiles on evaluation data
qrF_prediction <-   predict(qrF_model,
                            newdata = data_eval[, !(names(data_train) %in% c("date", "station", "obs"))],
                            what = qt_levels,
                            all = TRUE)

# Evaluation
ind_noNAinRow <- which(!apply(qrF_prediction, 1, anyNA))
ind_use <- intersect(which(!is.na(data_eval$obs)), ind_noNAinRow)

qrf_crps <- crps_sample(y = data_eval$obs[ind_use],
                        dat = qrF_prediction[ind_use,])

summary(qrf_crps)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1706  0.4732  0.6653  0.9492  1.1198 18.3247


# Variable importance
qrF_model$importance
importance(qrF_model, scale = TRUE)

plot(importance(qrF_model))

# t2m_mean       2768911.7878
# t2m_var          31306.5937
# cape_mean        58513.2436
# cape_var        107881.1488
# sp_mean          20301.6519
# sp_var           17134.2683
# tcc_mean         20192.2325
# tcc_var          13246.7371
# u_pl850_mean     31900.9411
# u_pl850_var      13689.0719
# v_pl850_mean     21149.5334
# v_pl850_var      13085.3433
# q_pl850_mean    949698.7449
# q_pl850_var      29890.9892
# u_pl500_mean     31220.8116
# u_pl500_var      14085.9995
# v_pl500_mean     20986.2375
# v_pl500_var      14705.8042
# gh_pl500_mean   147350.9623
# gh_pl500_var     20090.1830
# sshf_mean        26762.8232
# sshf_var         44440.7769
# slhf_mean       278381.2256
# slhf_var        305553.3006
# u10_mean         18747.9700
# u10_var          14354.2898
# v10_mean         20794.3926
# v10_var          13297.7998
# ssr_mean        169644.8408
# ssr_var          22463.9642
# str_mean         18781.2784
# str_var          12550.6565
# d2m_mean       1735668.6967
# d2m_var          15032.2295
# station_alt      89631.3051
# station_orog     20085.4631
# station_lsm        215.8518