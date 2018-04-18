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
pars_gl <- rbind(pars_gl, data.frame("m" = 51:100, "crps" = NA))

## local

m_val <- 1:100
fnames <-paste0(Rdata_dir, "alldata_climWindow_emos_local_m", m_val, ".Rdata")

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

pars_lo <- pars


# compare

plot_dir <- "/home/sebastian/Projects/PP_NN/code/standard_postprocessing/tuning_par_plots/"
pdf(paste0(plot_dir, "emos_alldata_climWindow_tuning.pdf"), width = 6, height = 5, pointsize = 12)
col_gl <- "darkorange"
col_lo <- "purple"
ylims <- range(na.omit(c(pars_gl$crps, pars_lo$crps)))
plot(crps ~ m, data = pars_gl, col = col_gl, type = "l",
     ylim = ylims)
lines(crps ~ m, data = pars_lo, col = col_lo)
# indicate nbest best choices
nbest <- 10
points(crps ~ m, data = pars_gl[order(pars_gl[,2])[1:nbest],],
       col = col_gl)
points(crps ~ m, data = pars_lo[order(pars_lo[,2])[1:nbest],],
       col = col_lo)
legend("topright", legend = c("global", "local"), 
       col = c(col_gl, col_lo), bty = "n", lty = 1)
dev.off()