## analysis code for boosting results / tuning parameter influence

## global

rm(list=ls())

Rdata_dir <- "/media/sebastian/Elements/Postproc_NN/model_data/boosting/"

maxit_try <- c(100, 250, 500, 1000, 2000, 5000)
nu_try <- c(0.1, 0.05, 0.2, 0.5)
mstop_try <- c("max", "aic", "bic", "cv")
link_try <-  c("log")

pars <- expand.grid(as.factor(maxit_try), 
                    as.factor(nu_try), 
                    as.factor(mstop_try), 
                    as.factor(link_try))
names(pars) <- c("maxit", "nu", "mstop", "link")

savenames <- paste0(Rdata_dir,
                    "boosting_global_maxit", pars$maxit,
                    "_nu", pars$nu,
                    "_mstop-", pars$mstop,
                    "_link-", pars$link,
                    ".Rdata")

res <- rep(NA, length(savenames))
for(file in savenames){
  if(file.exists(file)){
    load(file)
    res[which(savenames == file)] <- mean(crps_boosting, na.rm = TRUE)
  } else{
    res[which(savenames == file)] <- NA
  }
}

pars$crps <- res
pars

## NOTE: only log as link function appears to work, function crashes otherwise

pars <- subset(pars, link == "log")

head(pars[with(pars, order(crps)), ])

library(ggplot2)
library(gridExtra)

ggplot(aes(y = crps, x = maxit), data = pars) + geom_point()

p11 <- ggplot(pars, aes(maxit, crps, colour = nu)) + geom_point()
p12 <- ggplot(pars, aes(maxit, crps, colour = mstop)) + geom_point()

p21 <- ggplot(pars, aes(nu, crps, colour = maxit)) + geom_point()
p22 <- ggplot(pars, aes(nu, crps, colour = mstop)) + geom_point()

p31 <- ggplot(pars, aes(mstop, crps, colour = maxit)) + geom_point()
p32 <- ggplot(pars, aes(mstop, crps, colour = nu)) + geom_point()

pdf("boosting_global_tuning.pdf", width = 2*5, height = 3*5, pointsize = 12)
grid.arrange(p11, p12, 
             p21, p22, 
             p31, p32, 
             nrow = 3, ncol = 2)
dev.off()

# mstop does not matter
# larger maxit is better, particularly for small learning rate, but does no matter if >= 1000
# small learning rate is better (for large enough maxit)


## local

rm(list=ls())

Rdata_dir <- "/media/sebastian/Elements/Postproc_NN/model_data/boosting/"

maxit_try <- c(100, 250, 500, 1000, 2000, 5000)
nu_try <- c(0.1, 0.05, 0.2, 0.5)
mstop_try <- c("max", "aic", "bic", "cv")
link_try <-  c("log")

pars <- expand.grid(as.factor(maxit_try), 
                    as.factor(nu_try), 
                    as.factor(mstop_try), 
                    as.factor(link_try))
names(pars) <- c("maxit", "nu", "mstop", "link")

savenames <- paste0(Rdata_dir,
                    "boosting_local_maxit", pars$maxit,
                    "_nu", pars$nu,
                    "_mstop-", pars$mstop,
                    "_link-", pars$link,
                    ".Rdata")

res <- rep(NA, length(savenames))
for(file in savenames){
  if(file.exists(file)){
    load(file)
    res[which(savenames == file)] <- mean(crps_boosting, na.rm = TRUE)
  } else{
    res[which(savenames == file)] <- NA
  }
}

pars$crps <- res
pars

## NOTE: some strangely high CRPS values, only consider those < 10

pars <- subset(pars, crps < 1.2)

head(pars[with(pars, order(crps)), ])

library(ggplot2)
library(gridExtra)

ggplot(aes(y = crps, x = maxit), data = pars) + geom_point()

p11 <- ggplot(pars, aes(maxit, crps, colour = nu)) + geom_point()
p12 <- ggplot(pars, aes(maxit, crps, colour = mstop)) + geom_point()

p21 <- ggplot(pars, aes(nu, crps, colour = maxit)) + geom_point()
p22 <- ggplot(pars, aes(nu, crps, colour = mstop)) + geom_point()

p31 <- ggplot(pars, aes(mstop, crps, colour = maxit)) + geom_point()
p32 <- ggplot(pars, aes(mstop, crps, colour = nu)) + geom_point()

pdf("boosting_local_tuning.pdf", width = 2*5, height = 3*5, pointsize = 12)
grid.arrange(p11, p12, 
             p21, p22, 
             p31, p32, 
             nrow = 3, ncol = 2)
dev.off()

## similar conclusions to global model, but mstop appears to matter more
## particularly using max appears to be benefitial
## further, smaller learning rates appear to be important, otherwise resulting CRPS can be large


## globalST

rm(list=ls())

Rdata_dir <- "/media/sebastian/Elements/Postproc_NN/model_data/boosting/"

maxit_try <- c(250, 500, 1000, 2000, 5000)
nu_try <- c(0.1, 0.05, 0.2)
mstop_try <- c("max", "aic", "bic", "cv")
link_try <-  c("log")

pars <- expand.grid(as.factor(maxit_try), 
                    as.factor(nu_try), 
                    as.factor(mstop_try), 
                    as.factor(link_try))
names(pars) <- c("maxit", "nu", "mstop", "link")

savenames <- paste0(Rdata_dir,
                    "boosting_globalST_maxit", pars$maxit,
                    "_nu", pars$nu,
                    "_mstop-", pars$mstop,
                    "_link-", pars$link,
                    ".Rdata")

res <- rep(NA, length(savenames))
for(file in savenames){
  if(file.exists(file)){
    load(file)
    res[which(savenames == file)] <- mean(crps_boosting, na.rm = TRUE)
  } else{
    res[which(savenames == file)] <- NA
  }
}

pars$crps <- res
pars

pars <- subset(pars, link == "log")

head(pars[with(pars, order(crps)), ])

library(ggplot2)
library(gridExtra)

ggplot(aes(y = crps, x = maxit), data = pars) + geom_point()

p11 <- ggplot(pars, aes(maxit, crps, colour = nu)) + geom_point()
p12 <- ggplot(pars, aes(maxit, crps, colour = mstop)) + geom_point()

p21 <- ggplot(pars, aes(nu, crps, colour = maxit)) + geom_point()
p22 <- ggplot(pars, aes(nu, crps, colour = mstop)) + geom_point()

p31 <- ggplot(pars, aes(mstop, crps, colour = maxit)) + geom_point()
p32 <- ggplot(pars, aes(mstop, crps, colour = nu)) + geom_point()

pdf("boosting_globalST_tuning.pdf", width = 2*5, height = 3*5, pointsize = 12)
grid.arrange(p11, p12,
             p21, p22,
             p31, p32,
             nrow = 3, ncol = 2)
dev.off()
