## comparison of EMOS local and EMOS global results

rm(list=ls())

load("/home/sebastian/Projects/PP_NN/code/standard_postprocessing/preliminary_results/scores_global.Rdata")
crpsout_global <- crpsout_n
load("/home/sebastian/Projects/PP_NN/code/standard_postprocessing/preliminary_results/scores_local.Rdata")
crpsout_local <- crpsout_n

rm(crpsout_n)

summary(crpsout_global)
#   Min.  1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5874  0.8945  1.0198  1.0654  1.1827  4.4714 
summary(crpsout_local)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.4960  0.7673  0.9008  0.9578  1.0756  5.2003      73 
summary(crpsout_global - crpsout_local)
#      Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -0.85386  0.04866  0.11477  0.10693  0.17393  0.82847       73 

produce_pdfs <- TRUE
pdf_folder <- "/home/sebastian/Projects/PP_NN/code/standard_postprocessing/preliminary_results/"

dates_fit <- seq(as.POSIXct("2008-01-01 00:00", tz = "UTC", origin = "1970-01-01 00:00"), 
                 as.POSIXct("2016-12-31 00:00", tz = "UTC", origin = "1970-01-01 00:00"), 
                 by="+1 day")

if(produce_pdfs){
  pdf(paste0(pdf_folder, "crpsdiff_emos_global-vs-local.pdf"), width = 6, height = 5, pointsize = 12)
}
plot(dates_fit, crpsout_global - crpsout_local, type = "l")
abline(h = 0, lty = 2)
if(produce_pdfs){
  dev.off()
}


## rolling mean
# based on http://www.cookbook-r.com/Manipulating_data/Calculating_a_moving_average/
movingAverage_ignoreNA <- function(x, n = 1) {
  before <- floor((n-1)/2)
  after <- ceiling((n-1)/2)
  
  # Track the sum and count of number of non-NA items
  s <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new <- c(rep(NA, i), x[1:(length(x)-i)])
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}


# plot 2: rolling 50-day mean (mean) CRPS values and differences (smoothed version of plot 1)
k <- 50
if(produce_pdfs){
  pdf(paste0(pdf_folder, "crps_rollingmean_emos_global-vs-local.pdf"), width = 12, height = 5, pointsize = 12)
}
par(mfrow = c(1,2))
plot(dates_fit, movingAverage_ignoreNA(crpsout_global, k), type = "l", ylim = c(0.75,1.65),
     main = "50 day rolling mean", ylab = "mean CRPS", col = "darkorange")
lines(dates_fit, movingAverage_ignoreNA(crpsout_local, k), col = "purple")
lines(dates_fit, movingAverage_ignoreNA(crpsout_ens, k), col = "gray", lty = 1)
legend("topright", legend = c("EMOS global", "EMOS local", "Ensemble"), lty = c(1,1, 1), col = c("darkorange", "purple", "gray"), bty = "n")

plot(dates_fit, movingAverage_ignoreNA(crpsout_global, k)-movingAverage_ignoreNA(crpsout_local, k), type = "l", 
     main = "50 day rolling mean", ylab = "mean CRPS difference: global - local")
abline(h = 0, lty = 2)
if(produce_pdfs){
  dev.off()
}

# plot 3: separate forecast cases by month and average over all years
months_fit <- format(dates_fit, "%m")
crpsout_ens_monthly <- NULL
crpsout_global_monthly <- NULL
crpsout_local_monthly <- NULL
for(i in 1:length(unique(months_fit))){
  ind_use <- which(!is.na(crpsout_local))
  crpsout_ens_monthly[i] <- mean(crpsout_ens[intersect(ind_use, which(months_fit == unique(months_fit)[i]))])
  crpsout_global_monthly[i] <- mean(crpsout_global[intersect(ind_use, which(months_fit == unique(months_fit)[i]))])
  crpsout_local_monthly[i] <- mean(crpsout_local[intersect(ind_use, which(months_fit == unique(months_fit)[i]))])
}

if(produce_pdfs){
  pdf(paste0(pdf_folder, "crps_monthlymean_emos_global-vs-local.pdf"), width = 6, height = 5, pointsize = 12)
}
plot(unique(format(dates_fit, "%m")), crpsout_ens_monthly, type = "o", ylim = range(c(crpsout_ens_monthly, crpsout_global_monthly, crpsout_local_monthly)),
     main = "monthly mean CRPS", ylab = "mean CRPS", xlab = "month", col = "gray")
lines(unique(format(dates_fit, "%m")), crpsout_global_monthly, col = "darkorange", type = "o")
lines(unique(format(dates_fit, "%m")), crpsout_local_monthly, col = "purple", type = "o")
legend("top", legend = c("Ensemble", "EMOS global", "EMOS local"), lty = c(1,1,1), col = c("gray", "darkorange", "purple"), bty = "n")
if(produce_pdfs){
  dev.off()
}