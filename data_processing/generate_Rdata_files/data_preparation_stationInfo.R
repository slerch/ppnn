## extract station information and save as Rdata file for globalST models

rm(list=ls())

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_all.Rdata"))

library(ncdf4)
nc <- nc_open(paste0(data_dir, "data_aux_geo_interpolated.nc"))
nc
station_metadata <- list()
station_metadata$station_id <- ncvar_get(nc, "station_id")
station_metadata$lsm <- ncvar_get(nc, "lsm")
station_metadata$orog <- ncvar_get(nc, "orog")
station_metadata$station_alt <- ncvar_get(nc, "station_alt")
nc_close(nc)

any(unique(data$station) != station_metadata$station_id)
station_alt <- rep(station_metadata$station_alt, length(unique(data$date)))
station_orog <- rep(station_metadata$orog, length(unique(data$date)))
station_lsm <- as.factor(rep(station_metadata$lsm, length(unique(data$date))))

save(station_alt, station_orog, station_lsm,
     file = paste0(data_dir, "data_all_stationInfo.Rdata"))

## ---------------------------------------- ##

## read lon/lat/alt for evaluation data set

rm(list=ls())

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_all.Rdata"))

library(ncdf4)
nc <- nc_open(paste0(data_dir, "data_aux_geo_interpolated.nc"))
nc
station_metadata <- list()
station_metadata$station_id <- ncvar_get(nc, "station_id")
station_metadata$lat <- ncvar_get(nc, "station_lat")
station_metadata$lon <- ncvar_get(nc, "station_lon")
station_metadata$lsm <- ncvar_get(nc, "lsm")
station_metadata$orog <- ncvar_get(nc, "orog")
station_metadata$station_alt <- ncvar_get(nc, "station_alt")
nc_close(nc)

station_info <- as.data.frame(station_metadata)

head(station_info)

save(station_info,
     file = paste0(data_dir, "data_eval_stationInfo.Rdata"))
