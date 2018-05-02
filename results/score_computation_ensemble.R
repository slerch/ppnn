rm(list=ls())

library(lubridate)

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_ensfc.Rdata"))

start_eval <- as.Date("2016-01-01 00:00", tz = "UTC")
end_eval <- as.Date("2016-12-31 00:00", tz = "UTC") 
data_eval_all <- subset(data_ens, date >= start_eval & date <= end_eval)
rm(data_ens)

ind_use <- which(!is.na(data_eval_all$obs))
data_eval <- data_eval_all[ind_use,]

library(scoringRules)

fc_matrix <- as.matrix(data_eval[,4:ncol(data_eval)])
crps_ens <- crps_sample(y = data_eval$obs, dat = fc_matrix)
vrh_ens <- sapply(seq_along(1:nrow(fc_matrix)),
                  function(i) rank(c(data_eval$obs[i], fc_matrix[i,]))[1])
ae_ens <- abs(data_eval$obs - apply(fc_matrix, 1, median))


## cut to noNA_cases_training15 later
load("CRPS_AE_PIT_all.Rdata")

df_res_ens <- data_eval[noNA_cases_training15, 1:3]
df_res_ens$crps <- crps_ens[noNA_cases_training15]
df_res_ens$vrh <- vrh_ens[noNA_cases_training15]
df_res_ens$ae <- ae_ens[noNA_cases_training15]

save(df_res_ens,
     file = "CRPS_AE_PIT_ens.Rdata")
