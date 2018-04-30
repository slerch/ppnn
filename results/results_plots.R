rm(list=ls())

load("CRPS_AE_PIT_all.Rdata")
load("CRPS_AE_PIT_ens.Rdata")
load("/media/sebastian/Elements/Postproc_NN/data/data_eval_stationInfo.Rdata")

ps_plot <- 2

## ---------- map plot of station locations ---------- ##

library(ggmap)
library(colorspace)

# v1: only black points

map <- get_googlemap(
  center = c(10, 51), zoom = 6, maptype = "terrain", scale = 1,
  style = 'feature:road|element:all|visibility:off&style=feature:administrative.locality|element:labels|visibility:off'
)
p <- ggmap(map)
p <- p + geom_point(data = station_info, aes(x = lon, y = lat), size=ps_plot) 
p <- p + scale_x_continuous(limits = range(station_info$lon)) + 
  scale_y_continuous(limits = range(station_info$lat) + c(0,-0.17)) 
p <- p + xlab("Longitude") + ylab("Latitude")
p

# v2: point color depends on altitude

mypal <- gray.colors(25, start = 0.1, end = 0.9)
# for colored figures:
# mypal <- diverge_hcl(n = 25, power = 0.75)
# or sequential_hcl(...)

val_range <- range(station_info$station_alt)

map <- get_googlemap(
  center = c(10, 51), zoom = 6, maptype = "terrain", scale = 1,
  style = 'feature:road|element:all|visibility:off&style=feature:administrative.locality|element:labels|visibility:off'
)
p <- ggmap(map)
p <- p + geom_point(data = station_info, aes(x = lon, y = lat, color = station_alt), size=ps_plot) 
p <- p + scale_x_continuous(limits = range(station_info$lon)) + 
  scale_y_continuous(limits = range(station_info$lat)+ c(0,-0.17)) 
p <- p + scale_color_gradientn(colors = mypal, limits = val_range, name = "Altitude")
p <- p + xlab("Longitude") + ylab("Latitude")
p

pdf("map.pdf", width = 5, height = 5, pointsize = 12)
p
dev.off()

## ---------- compute station-specific results ---------- ##

df_stat_15 <- df_stat_0715 <- station_info[,c(1:3,6)] 
names(df_stat_15)[4] <- names(df_stat_0715)[4] <- "alt"

df_stat_15[,5:13] <- NA
for(j in 5:13){
  names(df_stat_15)[j] <- names(df_crps_15)[j-1]
}
df_stat_15[,14] <- NA
names(df_stat_15)[14] <- "ens"

df_stat_0715[,5:14] <- NA
for(j in 5:14){
  names(df_stat_0715)[j] <- names(df_crps_0715)[j-1]
}
df_stat_0715[,15] <- NA
names(df_stat_0715)[15] <- "ens"


for(i in 1:length(df_stat_15$station_id)){
  st <- df_stat_15$station_id[i]
  # 2015 pp models
  df_st_pp_15 <- subset(df_crps_15, station == st)
  if(nrow(df_st_pp_15) == 0){
    next
  }
  df_st_pp_15_mean <- apply(df_st_pp_15[,4:ncol(df_st_pp_15)], 2, mean)
  df_stat_15[i,5:13] <- df_st_pp_15_mean
  
  # 2007-2015 pp models
  df_st_pp_0715 <- subset(df_crps_0715, station == st)
  df_st_pp_0715_mean <- apply(df_st_pp_0715[,4:ncol(df_st_pp_0715)], 2, mean)
  df_stat_0715[i,5:14] <- df_st_pp_0715_mean
  
  # Ensemble
  df_st_ens <- subset(df_res_ens, station == st)
  df_stat_15[i,14] <- df_stat_0715[i,15] <- mean(df_st_ens$crps)
}


df_stat_15_noNA <- df_stat_15[complete.cases(df_stat_15),]
df_stat_0715_noNA <- df_stat_0715[complete.cases(df_stat_0715),]


## ---------- boxplots of station-specific mean CRPSS ---------- ##

compute_ss <- function(fc, ref){
  1-fc/ref
}

# relative to ensemble
df_crpss_ens_15 <- df_stat_15_noNA[,1:(ncol(df_stat_15_noNA)-1)]
fcnames <- names(df_stat_15_noNA[,5:(ncol(df_stat_15_noNA)-1)])
for(i in 1:length(fcnames)){
  df_crpss_ens_15[,4+i] <- compute_ss(df_stat_15_noNA[,fcnames[i]], df_stat_15_noNA$ens)
}

df_crpss_ens_0715 <- df_stat_0715_noNA[,1:(ncol(df_stat_0715_noNA)-1)]
fcnames <- names(df_stat_0715_noNA[,5:(ncol(df_stat_0715_noNA)-1)])
for(i in 1:length(fcnames)){
  df_crpss_ens_0715[,4+i] <- compute_ss(df_stat_0715_noNA[,fcnames[i]], df_stat_0715_noNA$ens)
}


ww <- 10
hh <- 6
ps <- 12

pdf("CRPSS_ens.pdf", width = ww, height = hh, pointsize = ps)
par(mfrow=c(1,2), 
    las = 2,
    mar = c(5,4,4,2) + 0.1 + c(2.5,0,0,0))

plotlim <- range(c(df_crpss_ens_15[,5:13]), df_crpss_ens_0715[,5:14])
plotlim <- c(-0.75, plotlim[2])
boxplot(df_crpss_ens_15[,5:13], main = "Training 2015", ylim = plotlim,
        pch = 20, cex = 0.5, axes = FALSE,
        ylab = "Mean CRPSS")
axis(2)
axis(1, at = 1:9, labels = c("FCN", "FCN-aux", "FCN-emb", "FCN-aux-emb", "NN-aux-emb", "EMOS-gl", "EMOS-loc", "EMOS-loc-boost", "QRF"))
# box()
abline(h=0, lty =2)

boxplot(df_crpss_ens_0715[,5:14], main = "Training 2007-2015", ylim = plotlim,
        pch = 20, cex = 0.5, axes = FALSE,
        ylab = "Mean CRPSS")
axis(2)
axis(1, at = 1:10, labels = c("FCN", "FCN-aux", "FCN-emb", "FCN-aux-emb", "NN-aux", "NN-aux-emb", "EMOS-gl", "EMOS-loc", "EMOS-loc-boost", "QRF"))
# box()
abline(h=0, lty =2)
dev.off()

# relative to EMOS-loc

df_crpss_emosloc_15 <- df_stat_15_noNA[,!(names(df_stat_15_noNA) == "emos_loc_15")]
fcnames <- names(df_crpss_emosloc_15)[5:ncol(df_crpss_emosloc_15)]
for(i in 1:length(fcnames)){
  df_crpss_emosloc_15[,4+i] <- compute_ss(df_stat_15_noNA[,fcnames[i]], df_stat_15_noNA$emos_loc_15)
}

df_crpss_emosloc_0715 <- df_stat_0715_noNA[,!(names(df_stat_0715_noNA) == "emos_loc_0715")]
fcnames <- names(df_crpss_emosloc_0715)[5:ncol(df_crpss_emosloc_0715)]
for(i in 1:length(fcnames)){
  df_crpss_emosloc_0715[,4+i] <- compute_ss(df_stat_0715_noNA[,fcnames[i]], df_stat_0715_noNA$emos_loc_0715)
}


# boxplot(df_crpss_emosloc_15[,5:12], main = "Training 2015")
# abline(h=0, lty =2)
# 
# par(mfrow=c(1,2))
# plotlim <- range(c(df_crpss_emosloc_15[,5:13]), df_crpss_emosloc_0715[,5:14])
# boxplot(df_crpss_emosloc_15[,5:13], main = "Training 2015", ylim = plotlim)
# abline(h=0, lty =2)
# boxplot(df_crpss_emosloc_0715[,5:14], main = "Training 2007-2015", ylim = plotlim)
# abline(h=0, lty =2)


pdf("CRPSS_emos-loc.pdf", width = ww, height = hh, pointsize = ps)

par(mfrow=c(1,2), 
    las = 2,
    mar = c(5,4,4,2) + 0.1 + c(2.5,0,0,0))

plotlim <- range(c(df_crpss_emosloc_15[,5:13]), df_crpss_emosloc_0715[,5:14])
plotlim <- c(-0.75, plotlim[2])
boxplot(df_crpss_emosloc_15[,5:13], main = "Training 2015", ylim = plotlim,
        pch = 20, cex = 0.5, axes = FALSE,
        ylab = "Mean CRPSS")
axis(2)
axis(1, at = 1:9, labels = c("FCN", "FCN-aux", "FCN-emb", "FCN-aux-emb", "NN-aux-emb", "EMOS-gl", "EMOS-loc-boost", "QRF", "Ensemble"))
# box()
abline(h=0, lty =2)

boxplot(df_crpss_emosloc_0715[,5:14], main = "Training 2007-2015", ylim = plotlim,
        pch = 20, cex = 0.5, axes = FALSE,
        ylab = "Mean CRPSS")
axis(2)
axis(1, at = 1:10, labels = c("FCN", "FCN-aux", "FCN-emb", "FCN-aux-emb", "NN-aux", "NN-aux-emb", "EMOS-gl", "EMOS-loc-boost", "QRF", "Ensemble"))
# box()
abline(h=0, lty =2)

dev.off()



# relative to boosting

df_crpss_bst_15 <- df_stat_15_noNA[,!(names(df_stat_15_noNA) == "bst_15")]
fcnames <- names(df_crpss_bst_15)[5:ncol(df_crpss_bst_15)]
for(i in 1:length(fcnames)){
  df_crpss_bst_15[,4+i] <- compute_ss(df_stat_15_noNA[,fcnames[i]], df_stat_15_noNA$bst_15)
}

df_crpss_bst_0715 <- df_stat_0715_noNA[,!(names(df_stat_0715_noNA) == "bst_0715")]
fcnames <- names(df_crpss_bst_0715)[5:ncol(df_crpss_bst_0715)]
for(i in 1:length(fcnames)){
  df_crpss_bst_0715[,4+i] <- compute_ss(df_stat_0715_noNA[,fcnames[i]], df_stat_0715_noNA$bst_0715)
}


# par(mfrow=c(1,2))
# plotlim <- range(c(df_crpss_bst_15[,5:13]), df_crpss_bst_0715[,5:14])
# plotlim <- c(-1,max(plotlim))
# 
# boxplot(df_crpss_bst_15[,5:13], main = "Training 2015", ylim = plotlim)
# abline(h=0, lty =2)
# boxplot(df_crpss_bst_0715[,5:14], main = "Training 2007-2015", ylim = plotlim)
# abline(h=0, lty =2)

pdf("CRPSS_emos-loc-boost.pdf", width = ww, height = hh, pointsize = ps)

par(mfrow=c(1,2), 
    las = 2,
    mar = c(5,4,4,2) + 0.1 + c(2.5,0,0,0))

plotlim <- range(c(df_crpss_bst_15[,5:13]), df_crpss_bst_0715[,5:14])
plotlim <- c(-1, plotlim[2])
boxplot(df_crpss_bst_15[,5:13], main = "Training 2015", ylim = plotlim,
        pch = 20, cex = 0.5, axes = FALSE,
        ylab = "Mean CRPSS")
axis(2)
axis(1, at = 1:9, labels = c("FCN", "FCN-aux", "FCN-emb", "FCN-aux-emb", "NN-aux-emb", "EMOS-gl", "EMOS-loc",  "QRF", "Ensemble"))
# box()
abline(h=0, lty =2)

boxplot(df_crpss_bst_0715[,5:14], main = "Training 2007-2015", ylim = plotlim,
        pch = 20, cex = 0.5, axes = FALSE,
        ylab = "Mean CRPSS")
axis(2)
axis(1, at = 1:10, labels = c("FCN", "FCN-aux", "FCN-emb", "FCN-aux-emb", "NN-aux", "NN-aux-emb", "EMOS-gl", "EMOS-loc", "QRF", "Ensemble"))
# box()
abline(h=0, lty =2)

dev.off()


# relative to qrf

df_crpss_qrf_15 <- df_stat_15_noNA[,!(names(df_stat_15_noNA) == "qrf_15")]
fcnames <- names(df_crpss_qrf_15)[5:ncol(df_crpss_qrf_15)]
for(i in 1:length(fcnames)){
  df_crpss_qrf_15[,4+i] <- compute_ss(df_stat_15_noNA[,fcnames[i]], df_stat_15_noNA$qrf_15)
}

df_crpss_qrf_0715 <- df_stat_0715_noNA[,!(names(df_stat_0715_noNA) == "qrf_0715")]
fcnames <- names(df_crpss_qrf_0715)[5:ncol(df_crpss_qrf_0715)]
for(i in 1:length(fcnames)){
  df_crpss_qrf_0715[,4+i] <- compute_ss(df_stat_0715_noNA[,fcnames[i]], df_stat_0715_noNA$qrf_0715)
}


# par(mfrow=c(1,2))
# plotlim <- range(c(df_crpss_qrf_15[,5:13]), df_crpss_qrf_0715[,5:14])
# plotlim <- c(-1,max(plotlim))
# 
# boxplot(df_crpss_qrf_15[,5:13], main = "Training 2015", ylim = plotlim)
# abline(h=0, lty =2)
# boxplot(df_crpss_qrf_0715[,5:14], main = "Training 2007-2015", ylim = plotlim)
# abline(h=0, lty =2)


pdf("CRPSS_qrf.pdf", width = ww, height = hh, pointsize = ps)

par(mfrow=c(1,2), 
    las = 2,
    mar = c(5,4,4,2) + 0.1 + c(2.5,0,0,0))

plotlim <- range(c(df_crpss_qrf_15[,5:13]), df_crpss_qrf_0715[,5:14])
plotlim <- c(-1, plotlim[2])
boxplot(df_crpss_qrf_15[,5:13], main = "Training 2015", ylim = plotlim,
        pch = 20, cex = 0.5, axes = FALSE,
        ylab = "Mean CRPSS")
axis(2)
axis(1, at = 1:9, labels = c("FCN", "FCN-aux", "FCN-emb", "FCN-aux-emb", "NN-aux-emb", "EMOS-gl", "EMOS-loc", "EMOS-loc-boost", "Ensemble"))
# box()
abline(h=0, lty =2)

boxplot(df_crpss_qrf_0715[,5:14], main = "Training 2007-2015", ylim = plotlim,
        pch = 20, cex = 0.5, axes = FALSE,
        ylab = "Mean CRPSS")
axis(2)
axis(1, at = 1:10, labels = c("FCN", "FCN-aux", "FCN-emb", "FCN-aux-emb", "NN-aux", "NN-aux-emb", "EMOS-gl", "EMOS-loc", "EMOS-loc-boost", "Ensemble"))
# box()
abline(h=0, lty =2)

dev.off()




# relative to nn_aux_emb

df_crpss_nn_aux_emb_15 <- df_stat_15_noNA[,!(names(df_stat_15_noNA) == "nn_aux_emb_15")]
fcnames <- names(df_crpss_nn_aux_emb_15)[5:ncol(df_crpss_nn_aux_emb_15)]
for(i in 1:length(fcnames)){
  df_crpss_nn_aux_emb_15[,4+i] <- compute_ss(df_stat_15_noNA[,fcnames[i]], df_stat_15_noNA$nn_aux_emb_15)
}

df_crpss_nn_aux_emb_0715 <- df_stat_0715_noNA[,!(names(df_stat_0715_noNA) == "nn_aux_emb_0715")]
fcnames <- names(df_crpss_nn_aux_emb_0715)[5:ncol(df_crpss_nn_aux_emb_0715)]
for(i in 1:length(fcnames)){
  df_crpss_nn_aux_emb_0715[,4+i] <- compute_ss(df_stat_0715_noNA[,fcnames[i]], df_stat_0715_noNA$nn_aux_emb_0715)
}


# par(mfrow=c(1,2))
# plotlim <- range(c(df_crpss_nn_aux_emb_15[,5:13]), df_crpss_nn_aux_emb_0715[,5:14])
# plotlim <- c(-1,max(plotlim))
# 
# boxplot(df_crpss_nn_aux_emb_15[,5:13], main = "Training 2015", ylim = plotlim)
# abline(h=0, lty =2)
# boxplot(df_crpss_nn_aux_emb_0715[,5:14], main = "Training 2007-2015", ylim = plotlim)
# abline(h=0, lty =2)


pdf("CRPSS_nn-aux-emb.pdf", width = ww, height = hh, pointsize = ps)

par(mfrow=c(1,2), 
    las = 2,
    mar = c(5,4,4,2) + 0.1 + c(2.5,0,0,0))

plotlim <- range(c(df_crpss_nn_aux_emb_15[,5:13]), df_crpss_nn_aux_emb_0715[,5:14])
plotlim <- c(-1, plotlim[2])
boxplot(df_crpss_nn_aux_emb_15[,5:13], main = "Training 2015", ylim = plotlim,
        pch = 20, cex = 0.5, axes = FALSE,
        ylab = "Mean CRPSS")
axis(2)
axis(1, at = 1:9, labels = c("FCN", "FCN-aux", "FCN-emb", "FCN-aux-emb",  "EMOS-gl", "EMOS-loc", "EMOS-loc-boost", "QRF", "Ensemble"))
# box()
abline(h=0, lty =2)

boxplot(df_crpss_nn_aux_emb_0715[,5:14], main = "Training 2007-2015", ylim = plotlim,
        pch = 20, cex = 0.5, axes = FALSE,
        ylab = "Mean CRPSS")
axis(2)
axis(1, at = 1:10, labels = c("FCN", "FCN-aux", "FCN-emb", "FCN-aux-emb", "NN-aux", "EMOS-gl", "EMOS-loc", "EMOS-loc-boost", "QRF", "Ensemble"))
# box()
abline(h=0, lty =2)

dev.off()




## ---------- plot of best model by station ---------- ##

df_best_15 <- df_best_0715 <- df_stat[,1:3]
df_best_15$model <- NA
df_best_0715$model <- NA
df_best_15$diff_to_2nd <- NA
df_best_0715$diff_to_2nd <- NA

for(i in 1:length(df_best_15$station_id)){
  st <- df_best_15$station_id[i]
  
  # 2015
  df_st <- subset(df_crps_15, station == st)
  if(nrow(df_st) == 0){
    next
  }
  df_st_means <- apply(df_st[,4:ncol(df_st)], 2, mean)
  best <- names(df_st_means)[which(df_st_means == sort(df_st_means)[1])]
  second <- names(df_st_means)[which(df_st_means == sort(df_st_means)[2])]
  diff <- df_st_means[best] - df_st_means[second]
  df_best_15$model[i] <- best
  df_best_15$diff_to_2nd[i] <- as.numeric(diff)
  
  # 2007-15
  df_st <- subset(df_crps_0715, station == st)
  df_st_means <- apply(df_st[,4:ncol(df_st)], 2, mean)
  best <- names(df_st_means)[which(df_st_means == sort(df_st_means)[1])]
  second <- names(df_st_means)[which(df_st_means == sort(df_st_means)[2])]
  diff <- df_st_means[best] - df_st_means[second]
  df_best_0715$model[i] <- best
  df_best_0715$diff_to_2nd[i] <- as.numeric(diff)
}


df_best_15_noNA <- df_best_15[complete.cases(df_best_15),]
df_best_0715_noNA <- df_best_0715[complete.cases(df_best_0715),]


library(colorspace)

allmodels_0715 <- names(df_crps_0715)[4:ncol(df_crps_0715)]
allmodels_15 <- names(df_crps_15)[4:ncol(df_crps_15)]
mypal <- rainbow_hcl(length(allmodels_0715))

mypal_use <- c("fc_15" = mypal[1],
               "fc_aux_15" = mypal[9],
               "fc_emb_15" = mypal[3],
               "fc_aux_emb_15" = mypal[4],
               "nn_aux_15" = mypal[5],
               "nn_aux_emb_15" = mypal[8],
               "emos_gl_15" = mypal[6],
               "emos_loc_15" = mypal[7],
               "bst_15" = mypal[2],
               "qrf_15" = mypal[10])

p <- ggmap(map)
p <- p + geom_point(data = df_best_15_noNA, aes(x = lon, y = lat, color = model), size=ps_plot)
p <- p + scale_x_continuous(limits = range(df_plot$lon)) +
  scale_y_continuous(limits = range(df_plot$lat)+ c(0,-0.17))
p <- p + xlab("Longitude") + ylab("Latitude")
p <- p + scale_colour_manual(values = mypal_use, name = "Best model",
                             labels = c("EMOS-loc-boost", "EMOS-gl", "EMOS-loc", "FCN-aux", "FCN-aux-emb", "FC-emb", "NN-aux", "NN-aux-emb", "QRF"),
                             limits = c("bst_15", "emos_gl_15", "emos_loc_15", "fc_aux_15", "fc_aux_emb_15",  "fc_emb_15", "nn_aux_15",   "nn_aux_emb_15", "qrf_15"))
p1 <- p

mypal_use <- c("fc_0715" = mypal[1],
               "fc_aux_0715" = mypal[9],
               "fc_emb_0715" = mypal[3],
               "fc_aux_emb_0715" = mypal[4],
               "nn_aux" = mypal[5],
               "nn_aux_emb_0715" = mypal[8],
               "emos_gl_0715" = mypal[6],
               "emos_loc_0715" = mypal[7],
               "bst_0715" = mypal[2],
               "qrf_0715" = mypal[10])

p <- ggmap(map)
p <- p + geom_point(data = df_best_0715_noNA, aes(x = lon, y = lat, color = model), size=ps_plot)
p <- p + scale_x_continuous(limits = range(df_plot$lon)) +
  scale_y_continuous(limits = range(df_plot$lat)+ c(0,-0.17))
p <- p + xlab("Longitude") + ylab("Latitude")
p <- p + scale_colour_manual(values = mypal_use, name = "Best model",
                             labels = c("EMOS-loc-boost", "EMOS-loc", "FCN-aux-emb", "NN-aux", "NN-aux-emb", "QRF"))
p2 <- p

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

legend <- g_legend(p1)

pdf("mapplot_bestmodel_v0.pdf", width = 10, height = 5, pointsize = 12)
grid.arrange(p1 + theme(legend.position = 'none', plot.title = element_text(hjust = 0.5)) + ggtitle("Training 2015"),
             p2 + theme(legend.position = 'none', plot.title = element_text(hjust = 0.5)) + ggtitle("Training 2007-2015"),
             legend,
             ncol=3, nrow=1, widths=c(3/7,3/7,1/7))
dev.off()

## alternative, with individual legends:

mypal <- rainbow_hcl(length(allmodels_0715))
mypal_use <- c("fc_15" = mypal[1], 
               "fc_aux_15" = mypal[9], 
               "fc_emb_15" = mypal[3], 
               "fc_aux_emb_15" = mypal[4], 
               "nn_aux_15" = mypal[5], 
               "nn_aux_emb_15" = mypal[8], 
               "emos_gl_15" = mypal[6], 
               "emos_loc_15" = mypal[7], 
               "bst_15" = mypal[2],
               "qrf_15" = mypal[10])

p <- ggmap(map)
p <- p + geom_point(data = df_best_15_noNA, aes(x = lon, y = lat, color = model), size=ps_plot)
p <- p + scale_x_continuous(limits = range(df_plot$lon)) +
  scale_y_continuous(limits = range(df_plot$lat)+ c(0,-0.17))
p <- p + xlab("Longitude") + ylab("Latitude")
p <- p + scale_colour_manual(values = mypal_use, name = "Best model",
                             labels = c("EMOS-loc-boost", "EMOS-gl", "EMOS-loc", "FCN-aux", "FCN-aux-emb", "FC-emb", "NN-aux-emb", "QRF"))
p <- p + ggtitle("Training 2015")  +
  theme(plot.title = element_text(hjust = 0.5))
p1 <- p

mypal_use <- c("fc_0715" = mypal[1], 
               "fc_aux_0715" = mypal[9], 
               "fc_emb_0715" = mypal[3], 
               "fc_aux_emb_0715" = mypal[4], 
               "nn_aux" = mypal[5], 
               "nn_aux_emb_0715" = mypal[8], 
               "emos_gl_0715" = mypal[6], 
               "emos_loc_0715" = mypal[7], 
               "bst_0715" = mypal[2],
               "qrf_0715" = mypal[10])

p <- ggmap(map)
p <- p + geom_point(data = df_best_0715_noNA, aes(x = lon, y = lat, color = model), size=ps_plot)
p <- p + scale_x_continuous(limits = range(df_plot$lon)) +
  scale_y_continuous(limits = range(df_plot$lat)+ c(0,-0.17))
p <- p + xlab("Longitude") + ylab("Latitude")
p <- p + scale_colour_manual(values = mypal_use, name = "Best model",
                             labels = c("EMOS-loc-boost", "EMOS-loc", "FCN-aux-emb", "NN-aux", "NN-aux-emb", "QRF"))
p <- p + ggtitle("Training 2007-2015")  +
  theme(plot.title = element_text(hjust = 0.5))
p2 <- p

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

legend1 <- g_legend(p1)
legend2 <- g_legend(p2)

pdf("mapplot_bestmodel_v1.pdf", width = 10, height = 5, pointsize = 12)
grid.arrange(p1 + theme(legend.position = 'none'),
             legend1,            
             p2 + theme(legend.position = 'none'),
             legend2,
             ncol=4, nrow=1, widths=c(2/8,1/8,2/8,1/8))
dev.off()


## ---------- PIT and VRH histograms ---------- ##

## 2015 training

str(df_pit_15)
df_pit_15$ens <- df_res_ens$vrh

pits <- list()
for(i in 1:10){
  
  pit_breaks <- seq(0, 1, 1/17)
  vrh_breaks <- seq(0.5, 51.5, 3)
  
  vec <- df_pit_15[,3+i]
  
  if(i <= 8){
    hh <- hist(vec, breaks = pit_breaks, plot = FALSE)
  }
  if(i > 8){
    hh <- hist(vec, breaks = vrh_breaks, plot = FALSE)
  }
  
  pits[[i]] <- hh
}

ymax <- max(pits[[10]]$density)
ncl <- 51
ymax_pit <- (ncl+1)*ymax

pdf("PIT_2015.pdf", width = 15, height = 20, pointsize = 12)
par(mfrow=c(4,3))
modelnames <- c("FCN", "FCN-aux", "FCN-emb", "FCN-aux-emb", "NN-aux-emb", "EMOS-gl", "EMOS-loc", "EMOS-loc-boost", "QRF", "Ensemble")
desired_order <- c("Ensemble", "EMOS-gl", "EMOS-loc", 
                   "EMOS-loc-boost", "QRF",
                   "FCN", "FCN-aux", "FCN-emb", 
                   "FCN-aux-emb", "NN-aux-emb")
for(i in 1:length(desired_order)){
  model <- desired_order[i]
  ind <- which(modelnames == model)
  if(ind > 8){
    plot(pits[[ind]], ylim = c(0,ymax), freq = FALSE, 
         col = "gray", border = "white",
         xlab = "Verification rank", main = model)
    abline(h = 1/51, lty = 2)
  }
  if(ind <= 8){
    plot(pits[[ind]], ylim = c(0,ymax_pit), freq = FALSE, 
         border = "white", col = "gray", 
         xlab = "PIT", main = model)
    abline(h = 1, lty = 2)
  }
  if(i == 5){
    plot.new() # start new row
  }
}
dev.off()


## 2007-2015 training

df_pit_0715$ens <- df_res_ens$vrh

pits <- list()
for(i in 1:11){
  
  pit_breaks <- seq(0, 1, 1/17)
  vrh_breaks <- seq(0.5, 51.5, 3)
  
  vec <- df_pit_0715[,3+i]
  
  if(i <= 9){
    hh <- hist(vec, breaks = pit_breaks, plot = FALSE)
  }
  if(i > 9){
    hh <- hist(vec, breaks = vrh_breaks, plot = FALSE)
  }
  
  pits[[i]] <- hh
}

ymax <- max(pits[[11]]$density)
ncl <- 51
ymax_pit <- (ncl+1)*ymax

pdf("PIT_2007-15.pdf", width = 15, height = 20, pointsize = 12)
par(mfrow=c(4,3))
modelnames <- c("FCN", "FCN-aux", "FCN-emb", "FCN-aux-emb", "NN-aux", "NN-aux-emb", "EMOS-gl", "EMOS-loc", "EMOS-loc-boost", "QRF", "Ensemble")
desired_order <- c("Ensemble", "EMOS-gl", "EMOS-loc", 
                   "EMOS-loc-boost", "QRF",
                   "FCN", "FCN-aux", "FCN-emb", 
                   "FCN-aux-emb", "NN-aux-emb", "NN-aux")
for(i in 1:length(desired_order)){
  model <- desired_order[i]
  ind <- which(modelnames == model)
  if(ind > 9){
    plot(pits[[ind]], ylim = c(0,ymax), freq = FALSE, 
         col = "gray", border = "white",
         xlab = "Verification rank", main = model)
    abline(h = 1/51, lty = 2)
  }
  if(ind <= 9){
    plot(pits[[ind]], ylim = c(0,ymax_pit), freq = FALSE, 
         border = "white", col = "gray", 
         xlab = "PIT", main = model)
    abline(h = 1, lty = 2)
  }
  if(i == 5){
    plot.new() # start new row
  }
}
dev.off()