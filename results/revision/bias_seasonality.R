rm(list=ls())

library(lubridate)

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_ensfc.Rdata"))

data_eval_all <- data_ens
ind_use <- which(!is.na(data_eval_all$obs))
data_eval <- data_eval_all[ind_use,]

library(scoringRules)

fc_matrix <- as.matrix(data_eval[,4:ncol(data_eval)])
bias_sample <- function(y, dat){
  apply(dat, 1, mean) - y
}
bias_ens <- bias_sample(y = data_eval$obs, dat = fc_matrix)

df <- data.frame(date = data_eval$date, bias = bias_ens)

df$year <- year(df$date)
df$month <- month(df$date)

bias_monthlymean <- matrix(NA, ncol = 12, nrow = length(c(2007:2016)))
bias_monthlyqlow <- matrix(NA, ncol = 12, nrow = length(c(2007:2016)))
bias_monthlyqhigh <- matrix(NA, ncol = 12, nrow = length(c(2007:2016)))

for(yy in c(2007:2016)){
  print(yy)
  dfy <- subset(df, year == yy)
  mean_mm <- NULL
  # low_mm <- NULL
  # high_mm <- NULL
  for(mm in 1:12){
    dfy_mm <- subset(dfy, month == mm)
    mean_mm[mm] <- mean(dfy_mm$bias)
    # low_mm[mm] <- quantile(dfy_mm$bias, 0.05)
    # high_mm[mm] <- quantile(dfy_mm$bias, 0.95)
  }
  bias_monthlymean[which(c(2007:2016) == yy),] <- mean_mm
  # bias_monthlyqlow[which(c(2007:2016) == yy),] <- low_mm
  # bias_monthlyqhigh[which(c(2007:2016) == yy),] <- high_mm
}


pdf("bias_seasonality.pdf", width = 10, height = 5, pointsize = 10)

par(mfrow=c(1,2))

colors <- rainbow(length(c(2007:2016)))
plot(bias_monthlymean[1,], type = "l", col = colors[1], ylim = range(c(bias_monthlymean)+c(0,0.5)),
     xlab = "Month", ylab = "Monthly mean bias of ensemble", axes = FALSE)
axis(1, at = 1:12, labels = month.abb)
axis(2)
for(row in 2:nrow(bias_monthlymean)){
  lines(bias_monthlymean[row,], col = colors[row])
}
legend("top", col = colors, lty = 1, legend = c(2007:2016), ncol = 4, bty = "n")



indpl <- 9:10
plot(bias_monthlymean[indpl,][1,], type = "l", col = colors[indpl][1], ylim = range(c(bias_monthlymean))+c(0,0.5),
     xlab = "Month", ylab = "Monthly mean bias of ensemble", axes = FALSE)
axis(1, at = 1:12, labels = month.abb)
axis(2)
for(row in 2:nrow(bias_monthlymean[indpl,])){
  lines(bias_monthlymean[indpl,][row,], col = colors[indpl][row])
}
legend("top", col = colors[indpl], lty = 1, legend = c(2007:2016)[indpl], ncol = 3, bty = "n")

dev.off()

