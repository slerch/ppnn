## analysis code for qrf results / tuning parameter influence

rm(list=ls())

Rdata_dir <- "/media/sebastian/Elements/Postproc_NN/model_data/qrf/cluster_runs_local/"

ntree.try <- c(125,250,500,1000,2500)
nodesize.try <- c(5,10,15,20,50)
mtry.try <- c(10,15,20,25)
replace.try <- c(0,1)

pars <- expand.grid(as.factor(ntree.try), 
                    as.factor(nodesize.try), 
                    as.factor(mtry.try), 
                    as.factor(replace.try))
names(pars) <- c("ntree", "nodesize", "mtry", "replace")

savenames <- paste0(Rdata_dir, "qrf_local_ntree", pars[,1], "_nodesize", pars[,2], "_mtry", pars[,3], "_repl", pars[,4], ".Rdata")

res <- rep(NA, length(savenames))
for(file in savenames){
  if(file.exists(file)){
    load(file)
    res[which(savenames == file)] <- mean(qrf_crps)
  }
}

pars$crps <- res
head(pars)

pars[with(pars, order(crps)), ]

## Plots

library(ggplot2)
library(gridExtra)
# plot(pars$crps)

# ggplot(aes(y = crps, x = ntree), data = pars) + geom_point()
p11 <- ggplot(pars, aes(ntree, crps, colour = nodesize)) + geom_point()
p12 <- ggplot(pars, aes(ntree, crps, colour = mtry)) + geom_point()
p13 <- ggplot(pars, aes(ntree, crps, colour = replace)) + geom_point()

# ggplot(aes(y = crps, x = nodesize), data = pars) + geom_point()
p21 <- ggplot(pars, aes(nodesize, crps, colour = ntree)) + geom_point()
p22 <- ggplot(pars, aes(nodesize, crps, colour = mtry)) + geom_point()
p23 <- ggplot(pars, aes(nodesize, crps, colour = replace)) + geom_point()

# ggplot(aes(y = crps, x = mtry), data = pars) + geom_point()
p31 <- ggplot(pars, aes(mtry, crps, colour = ntree)) + geom_point()
p32 <- ggplot(pars, aes(mtry, crps, colour = nodesize)) + geom_point()
p33 <- ggplot(pars, aes(mtry, crps, colour = replace)) + geom_point()

# ggplot(aes(y = crps, x = replace), data = pars) + geom_point()
# ggplot(pars, aes(replace, crps, colour = ntree)) + geom_point()
# ggplot(pars, aes(replace, crps, colour = nodesize)) + geom_point()
# ggplot(pars, aes(replace, crps, colour = mtry)) + geom_point()

pdf("qrf_local_tuning.pdf", width = 14, height = 14, pointsize = 12)
grid.arrange(p11, p12, p13,
             p21, p22, p23,
             p31, p32, p33,
             nrow = 3, ncol = 3)
dev.off()

## restrict attention to replace == 1 option

pars1 <- subset(pars, replace == 1)
ggplot(pars1, aes(ntree, crps, colour = nodesize)) + geom_point()
ggplot(pars1, aes(ntree, crps, colour = mtry)) + geom_point()

p11 <- ggplot(pars1, aes(ntree, crps, colour = nodesize)) + geom_point()
p12 <- ggplot(pars1, aes(ntree, crps, colour = mtry)) + geom_point()

p22 <- ggplot(pars1, aes(nodesize, crps, colour = mtry)) + geom_point()
p23 <- ggplot(pars1, aes(nodesize, crps, colour = replace)) + geom_point()

p31 <- ggplot(pars1, aes(mtry, crps, colour = ntree)) + geom_point()
p32 <- ggplot(pars1, aes(mtry, crps, colour = nodesize)) + geom_point()


pdf("qrf_local_tuning_replace1.pdf", width = 10, height = 15, pointsize = 12)
grid.arrange(p11, p12, 
             p21, p22, 
             p31, p32, 
             nrow = 3, ncol = 2)
dev.off()


# some impressions from plots:
# larger ntree benefital, but only small effect
# smaller nodesize is better, stronger effect
# larger mtry is better, also larger effect
# better with replacement, but effect not that large