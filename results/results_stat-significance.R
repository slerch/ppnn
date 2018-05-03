rm(list=ls())

load("CRPS_AE_PIT_all.Rdata")
load("CRPS_AE_PIT_ens.Rdata")
load("/media/sebastian/Elements/Postproc_NN/data/data_eval_stationInfo.Rdata")
nstations <- length(station_info$station_id)

## assess statistical significance of differences on a station level
## using Diebold-Mariano tests
## and applying a Benjamini-Hochberg-correction to p-values 
## (as suggested by Wilks 2016, BAMS)

## ---------- 2015 training ---------- ##

# compute station-wise test statistics for all model comparisons

df_crps_15$ens <- df_res_ens$crps
models_15 <- names(df_crps_15)[4:ncol(df_crps_15)]
nmod <- length(models_15)

# for each station, there is a nmod x nomod matrix of test statistic values to be saved
# dimensions = station, model1, model2
dm_tn <- array(NA, dim = c(nstations, nmod, nmod))

library(forecast)

for(i in 1:nstations){
  print(i)
  st <- station_info$station_id[i]
  
  # 2015 pp models
  dfst_15 <- subset(df_crps_15, station == st)
  if(nrow(dfst_15) == 0){
    next
  }
  
  for(m1 in 1:nmod){
    model1 <- models_15[m1]
    fc1 <- dfst_15[,model1]
    
    for(m2 in 1:nmod){
      model2 <- models_15[m2]
      fc2 <- dfst_15[,model2]
      
      if(!any(fc1 != fc2)){
        next
      }
      tmp <- dm.test(e1 = fc1, e2 = fc2,
                     h = 2, power = 1)
      dm_tn[i,m1,m2] <- tmp$statistic 
    }
  }
}

omit <- which(is.na(dm_tn[,1,2]))
dm_tn_withNA <- dm_tn
dm_tn <- dm_tn_withNA[-omit,,] 

freq_raw <- matrix(NA, nrow = nmod, ncol = nmod)
freq_BH <- matrix(NA, nrow = nmod, ncol = nmod)

bh_correction <- function(pval, level = 0.05){
  nn <- length(pval)
  pv_srt <- sort(pval)
  limits <- 1:nn*level/nn
  if(length(which(pv_srt < limits)) == 0){
    return(discoveries <- 0)
  } else{
    discoveries <- max(which(pv_srt < limits)) 
  }
  return(discoveries/nn)
}

level <- 0.05
for(m1 in 1:nmod){
  for(m2 in 1:nmod){
    if(m1 == m2){
      next
    } else if(m2 > m1){
      tn <- dm_tn[,m1,m2]
      pval <- pnorm(tn)
    } else if(m1 > m2){
      tn <- dm_tn[,m2,m1]
      pval <- 1-pnorm(tn)
    }
    freq_raw[m1,m2] <- sum(pval <= level)/length(pval)
    freq_BH[m1,m2] <- bh_correction(pval)
  }
}

freq_raw_15 <- freq_raw
freq_BH_15 <- freq_BH 
df_BH_15 <- data.frame(100*freq_BH_15)
names(df_BH_15) <- models_15
rownames(df_BH_15) <- models_15
df_BH_15

desired_order <- c("ens", "emos_gl_15", "emos_loc_15", "bst_15", "qrf_15",
                   "fc_15", "fc_aux_15", "fc_emb_15", "fc_aux_emb_15", "nn_aux_emb_15")
df_BH_15[desired_order,desired_order]

library(xtable)
xtable(df_BH_15[desired_order,desired_order], digits = 1)

## ---------- 2007-2015 training ---------- ##

# compute station-wise test statistics for all model comparisons

df_crps_0715$ens <- df_res_ens$crps
models_0715 <- names(df_crps_0715)[4:ncol(df_crps_0715)]
nmod <- length(models_0715)

# for each station, there is a nmod x nomod matrix of test statistic values to be saved
# dimensions = station, model1, model2
dm_tn <- array(NA, dim = c(nstations, nmod, nmod))

library(forecast)

for(i in 1:nstations){
  print(i)
  st <- station_info$station_id[i]
  
  # 2015 pp models
  dfst <- subset(df_crps_0715, station == st)
  if(nrow(dfst) == 0){
    next
  }
  
  for(m1 in 1:nmod){
    model1 <- models_0715[m1]
    fc1 <- dfst[,model1]
    
    for(m2 in 1:nmod){
      model2 <- models_0715[m2]
      fc2 <- dfst[,model2]
      
      if(!any(fc1 != fc2)){
        next
      }
      tmp <- dm.test(e1 = fc1, e2 = fc2,
                     h = 2, power = 1)
      dm_tn[i,m1,m2] <- tmp$statistic 
    }
  }
}

omit <- which(is.na(dm_tn[,1,2]))
dm_tn_withNA <- dm_tn
dm_tn <- dm_tn_withNA[-omit,,] 

freq_raw <- matrix(NA, nrow = nmod, ncol = nmod)
freq_BH <- matrix(NA, nrow = nmod, ncol = nmod)

bh_correction <- function(pval, level = 0.05){
  nn <- length(pval)
  pv_srt <- sort(pval)
  limits <- 1:nn*level/nn
  if(length(which(pv_srt < limits)) == 0){
    return(discoveries <- 0)
  } else{
    discoveries <- max(which(pv_srt < limits)) 
  }
  return(discoveries/nn)
}

level <- 0.05
for(m1 in 1:nmod){
  for(m2 in 1:nmod){
    if(m1 == m2){
      next
    } else if(m2 > m1){
      tn <- dm_tn[,m1,m2]
      pval <- pnorm(tn)
    } else if(m1 > m2){
      tn <- dm_tn[,m2,m1]
      pval <- 1-pnorm(tn)
    }
    freq_raw[m1,m2] <- sum(pval <= level)/length(pval)
    freq_BH[m1,m2] <- bh_correction(pval)
  }
}

freq_raw_0715 <- freq_raw
freq_BH_0715 <- freq_BH 
df_BH_0715 <- data.frame(100*freq_BH_0715)
names(df_BH_0715) <- models_0715
rownames(df_BH_0715) <- models_0715
df_BH_0715

desired_order <- c("ens", "emos_gl_0715", "emos_loc_0715", "bst_0715", "qrf_0715",
                   "fc_0715", "fc_aux_0715", "fc_emb_0715", "fc_aux_emb_0715", "nn_aux", "nn_aux_emb_0715")
df_BH_0715[desired_order,desired_order]

library(xtable)
xtable(df_BH_0715[desired_order,desired_order], digits = 1)
