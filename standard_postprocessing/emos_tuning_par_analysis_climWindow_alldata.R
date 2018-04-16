## global

rm(list=ls())

Rdata_dir <- "/media/sebastian/Elements/Postproc_NN/model_data/emos/"

m_val <- m_val <- 1:50
fnames <-paste0(Rdata_dir, "alldata_climWindow_emos_global_m", m_val, ".Rdata")

crps_out <- NULL
for(ff in fnames){
  if(file.exists(ff)){
    load(ff)
    crps_out[which(fnames == ff)] <- mean(crps_pp, na.rm = T)
  } else{
    crps_out[which(fnames == ff)] <- NA
  }
}

pars <- data.frame("m" = m_val, "crps" = crps_out)

plot(crps ~ m, data = pars, type = "l")

pars[order(pars[,2]),] # 22, 21, 20, ...

pars_gl <- pars

## local

# running.