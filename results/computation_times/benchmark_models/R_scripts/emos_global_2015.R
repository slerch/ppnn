cat(paste("Now running: EMOS global (training 2007-2015), started", Sys.time(), "\n"))
cat(paste("... should take less than 1 min", "\n"))

mydata <- data[, -which(!(names(data) %in% c("obs", "date", "station", "t2m_mean", "t2m_var")))]

objective_fun <- function(par, ensmean, ensvar, obs){
  m <- cbind(1, ensmean) %*% par[1:2]
  ssq_tmp <- cbind(1, ensvar) %*% par[3:4]
  if(any(ssq_tmp < 0)){
    return(999999)
  } else{
    s <- sqrt(ssq_tmp)
    return(sum(crps_norm(y = obs, location = m, scale = s)))
  }
}

gradfun_wrapper <- function(par, obs, ensmean, ensvar){
  loc <- cbind(1, ensmean) %*% par[1:2]
  sc <- sqrt(cbind(1, ensvar) %*% par[3:4])
  dcrps_dtheta <- gradcrps_norm(y = obs, location = loc, scale = sc) 
  out1 <- dcrps_dtheta[,1] %*% cbind(1, ensmean)
  out2 <- dcrps_dtheta[,2] %*% 
    cbind(1/(2*sqrt(par[3]+par[4]*ensvar)), 
          ensvar/(2*sqrt(par[3]+par[4]*ensvar)))
  return(as.numeric(cbind(out1,out2)))
}

train_end <- as.Date("2016-01-01 00:00 UTC") - days(2)
train_start <- as.Date("2015-01-01 00:00 UTC") 

data_train <- subset(mydata, date >= train_start & date <= train_end)
data_train <- data_train[complete.cases(data_train), ]


eval_start <- as.Date("2016-01-01 00:00 UTC")
eval_end <- as.Date("2016-12-31 00:00 UTC")
eval_dates <- seq(eval_start, eval_end, by = "1 day")

data_eval_all <- subset(mydata, date >= eval_start & date <= eval_end)

t1 <- Sys.time()
optim_out <- optim(par = c(1,1,1,1), 
                   fn = objective_fun,
                   gr = gradfun_wrapper,
                   ensmean = data_train$t2m_mean, 
                   ensvar = data_train$t2m_var, 
                   obs = data_train$obs,
                   method = "BFGS")

par_out <- optim_out$par

out_loc <- NULL
out_sc <- NULL

out_loc <- c(cbind(1, data_eval_all$t2m_mean) %*% par_out[1:2])
scsquared_tmp <- c(cbind(1, data_eval_all$t2m_var) %*% par_out[3:4])
out_sc <- sqrt(abs(scsquared_tmp))


t2 <- Sys.time()
time_elapsed <- difftime(t2, t1, units='mins')
cat("... finished, took", round(time_elapsed,2), "minutes", "\n")