rm(list=ls())

require(scoringRules)
if(packageVersion("scoringRules") < "0.9.4"){
  print("scoringRules package version 0.9.4 or newer is required")
}
require(lubridate)
require(quantregForest)
require(crch)

## Set parent directory
# should contain data_all.Rdata, and folder R_scripts/
dir <- "/media/sebastian/Elements/tmp/"
load(paste0(dir, "data_all.Rdata"))
data$sm_mean <- NULL
data$sm_var <- NULL

## training 2015 models

# EMOS global
source(paste0(dir, "R_scripts/emos_global_2015.R"))
time_emos_global_15 <- time_elapsed
rm(list=setdiff(ls(), c("dir", "data", 
                        "time_emos_global_15")))

# EMOS local, will take < 1 minute
source(paste0(dir, "R_scripts/emos_local_2015.R"))
time_emos_local_15 <- time_elapsed
rm(list=setdiff(ls(), c("dir", "data", 
                        "time_emos_global_15", "time_emos_local_15")))

# Boosting, will take around 15 minutes
source(paste0(dir, "R_scripts/boosting_2015.R"))
rm(list=setdiff(ls(), c("dir", "data", 
                        "time_emos_global_15", "time_emos_local_15", 
                        "time_bst_15")))

# QRF, will take around 8 minutes
source(paste0(dir, "R_scripts/qrf_2015.R"))
time_qrf_15 <- time_elapsed
rm(list=setdiff(ls(), c("dir", "data", 
                        "time_emos_global_15", "time_emos_local_15", 
                        "time_bst_15", "time_qrf_15")))

times_15 <- data.frame("emos_global_15" = time_emos_global_15, 
                       "emos_local_15" = time_emos_local_15,
                       "boosting_15" = time_bst_15,
                       "qrf_15" = time_qrf_15)

## 2007-15 Training

# EMOS global, will take < 1 minute
source(paste0(dir, "R_scripts/emos_global_2007-15.R"))
time_emos_global_0715 <- time_elapsed
rm(list=setdiff(ls(), c("dir", "data", 
                        "time_emos_global_15", "time_emos_local_15", 
                        "time_bst_15", "time_qrf_15",
                        "time_emos_global_0715")))

# EMOS local, will take < 1 minute
source(paste0(dir, "R_scripts/emos_local_2007-15.R"))
time_emos_local_0715 <- time_elapsed
rm(list=setdiff(ls(), c("dir", "data", 
                        "time_emos_global_15", "time_emos_local_15", 
                        "time_bst_15", "time_qrf_15",
                        "time_emos_global_0715", "time_emos_local_0715")))

# Boosting, will take around 45 minutes
source(paste0(dir, "R_scripts/boosting_2007-15.R"))
time_boosting_0715 <- time_elapsed
rm(list=setdiff(ls(), c("dir", "data", 
                        "time_emos_global_15", "time_emos_local_15", 
                        "time_bst_15", "time_qrf_15",
                        "time_emos_global_0715", "time_emos_local_0715",
                        "time_boosting_0715")))

# QRF, will take around 7 hours (!)
source(paste0(dir, "R_scripts/qrf_2007-15.R"))
time_qrf_0715 <- time_elapsed
rm(list=setdiff(ls(), c("dir", "data", 
                        "time_emos_global_15", "time_emos_local_15", 
                        "time_bst_15", "time_qrf_15",
                        "time_emos_global_0715", "time_emos_local_0715",
                        "time_boosting_0715", "time_qrf_0715")))

times_0715 <- data.frame("emos_global_0715" = time_emos_global_0715, 
                         "emos_local_0715" = time_emos_local_0715,
                         "boosting_0715" = time_bst_0715,
                         "qrf_0715" = time_qrf_0715)

##
## Results
##

## Training 2015

cat("Overview, training 2015", "\n")
times_15

## Training 2007-2015
cat("Overview, training 2007-2015", "\n")
times_0715