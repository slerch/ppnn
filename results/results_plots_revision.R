## plot of best model by station for revision of paper
## changes: remove country names, set background to B&W

rm(list=ls())

load("/media/sebastian/Elements/Postproc_NN/CRPS_AE_PIT_all.Rdata")
load("/media/sebastian/Elements/Postproc_NN/CRPS_AE_PIT_ens.Rdata")
load("/media/sebastian/Elements/Postproc_NN/data/data_eval_stationInfo.Rdata")

ps_plot <- 2

## ---------- map plot of station locations ---------- ##

library(ggmap)
library(colorspace)

## ---------- compute station-specific results ---------- ##

df_stat_15 <- df_stat_0715 <- station_info[,c(1:3,6)] 
names(df_stat_15)[4] <- names(df_stat_0715)[4] <- "alt"

df_stat_15[,5:14] <- NA
for(j in 5:14){
  names(df_stat_15)[j] <- names(df_crps_15)[j-1]
}
df_stat_15[,15] <- NA
names(df_stat_15)[15] <- "ens"

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
  df_stat_15[i,5:14] <- df_st_pp_15_mean
  
  # 2007-2015 pp models
  df_st_pp_0715 <- subset(df_crps_0715, station == st)
  df_st_pp_0715_mean <- apply(df_st_pp_0715[,4:ncol(df_st_pp_0715)], 2, mean)
  df_stat_0715[i,5:14] <- df_st_pp_0715_mean
  
  # Ensemble
  df_st_ens <- subset(df_res_ens, station == st)
  df_stat_15[i,15] <- df_stat_0715[i,15] <- mean(df_st_ens$crps)
}


df_stat_15_noNA <- df_stat_15[complete.cases(df_stat_15),]
df_stat_0715_noNA <- df_stat_0715[complete.cases(df_stat_0715),]

## ----- plot code

df_best_15 <- df_best_0715 <- df_stat_15[,1:3]
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

myshapes_use <- c("fc_15" = 15,
                  "fc_aux_15" = 15,
                  "fc_emb_15" = 15,
                  "fc_aux_emb_15" = 15,
                  "nn_aux_15" = 19,
                  "nn_aux_emb_15" = 19,
                  "emos_gl_15" = 2,
                  "emos_loc_15" = 2,
                  "bst_15" = 17,
                  "qrf_15" = 4) # 18

mylabels <- c("EMOS-gl", "EMOS-loc", "EMOS-loc-bst", "QRF" ,"FCN-aux", "FCN-emb", "FCN-aux-emb",  "NN-aux", "NN-aux-emb")
mybreaks <- c("emos_gl_15", "emos_loc_15", "bst_15", "qrf_15", "fc_aux_15", "fc_emb_15", "fc_aux_emb_15",   "nn_aux_15",   "nn_aux_emb_15")

map <- get_googlemap(
  center = c(10, 51), zoom = 6, maptype = "terrain", scale = 1,
  style = 'feature:road|element:all|visibility:off&style=feature:all|element:labels|visibility:off&sensor=false',
  color = "bw"
)

p <- ggmap(map)
p <- p + geom_point(data = df_best_15_noNA, 
                    aes(x = lon, y = lat, color = model, shape = model), 
                    size=ps_plot)
p <- p + scale_shape_manual(values=myshapes_use,
                            labels = mylabels,
                            limits = mybreaks,
                            name = "Best model")
p <- p + scale_colour_manual(values = mypal_use, 
                             labels = mylabels,
                             limits = mybreaks,
                             name = "Best model")
p <- p + scale_x_continuous(limits = range(df_best_15_noNA$lon)) +
  scale_y_continuous(limits = range(df_best_15_noNA$lat)+ c(0,-0.17))
p <- p + xlab("Longitude") + ylab("Latitude")
p1 <- p 


mypal_use <- c("fc_0715" = mypal[1],
               "fc_aux_0715" = mypal[9],
               "fc_emb_0715" = mypal[3],
               "fc_aux_emb_0715" = mypal[4],
               "nn_aux_0715" = mypal[5],
               "nn_aux_emb_0715" = mypal[8],
               "emos_gl_0715" = mypal[6],
               "emos_loc_0715" = mypal[7],
               "bst_0715" = mypal[2],
               "qrf_0715" = mypal[10])


myshapes_use <- c("fc_0715" = 15,
                  "fc_aux_0715" = 15,
                  "fc_emb_0715" = 15,
                  "fc_aux_emb_0715" = 15,
                  "nn_aux_0715" = 19,
                  "nn_aux_emb_0715" = 19,
                  "emos_gl_0715" = 2,
                  "emos_loc_0715" = 2,
                  "bst_0715" = 17,
                  "qrf_0715" = 4) # 18

mylabels <- c("EMOS-gl", "EMOS-loc", "EMOS-loc-bst", "QRF" ,"FCN-aux", "FCN-emb", "FCN-aux-emb","NN-aux", "NN-aux-emb")
mybreaks <- c("emos_gl_0715", "emos_loc_0715", "bst_0715", "qrf_0715", "fc_aux_0715", "fc_emb_0715", "fc_aux_emb_0715",  "nn_aux_0715",   "nn_aux_emb_0715")

p <- ggmap(map)
p <- p + geom_point(data = df_best_0715_noNA, 
                    aes(x = lon, y = lat, color = model, shape = model), 
                    size=ps_plot)
p <- p + scale_shape_manual(values=myshapes_use,
                            labels = mylabels,
                            limits = mybreaks,
                            name = "Best model")
p <- p + scale_colour_manual(values = mypal_use, 
                             labels = mylabels,
                             limits = mybreaks,
                             name = "Best model")
p <- p + scale_x_continuous(limits = range(df_best_0715_noNA$lon)) +
  scale_y_continuous(limits = range(df_best_0715_noNA$lat)+ c(0,-0.17))
p <- p + xlab("Longitude") + ylab("Latitude")
p2 <- p

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

legend <- g_legend(p1)
library(gridExtra)

pdf("mapplot_bestmodel_rev.pdf", width = 10, height = 6, pointsize = 12)
grid.arrange(p1 + theme(legend.position = 'none', plot.title = element_text(hjust = 0.5)) + ggtitle("Training 2015"),
             p2 + theme(legend.position = 'none', plot.title = element_text(hjust = 0.5)) + ggtitle("Training 2007-2015"),
             legend,
             ncol=3, nrow=1, widths=c(3,3,1))
dev.off()

## code for generating p1 may need to be run several times, otherwise all points are missing