## analysis code for qrf results / tuning parameter influence

rm(list=ls())

Rdata_dir <- "/media/sebastian/Elements/Postproc_NN/model_data/qrf/cluster_runs_global/"

# ntree.try <- c(125,250,500,1000,2500)
# nodesize.try <- c(5,10,15,20,50)
# mtry.try <- c(10,15,20,25)
# replace.try <- c(0,1)

# memTest version [+ pre-existing files]
ntree.try <- c(50, 100, 125, 250)
nodesize.try <- c(5, 10, 20, 50)
mtry.try <- c(10, 25)
replace.try <- c(1)

pars <- expand.grid(as.factor(ntree.try), 
                    as.factor(nodesize.try), 
                    as.factor(mtry.try), 
                    as.factor(replace.try))
names(pars) <- c("ntree", "nodesize", "mtry", "replace")

savenames <- paste0(Rdata_dir, "qrf_global_ntree", pars[,1], "_nodesize", pars[,2], "_mtry", pars[,3], "_repl", pars[,4], ".Rdata")

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
