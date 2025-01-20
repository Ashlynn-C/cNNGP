# set working directory to source file location
file.dir = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(file.dir))
rstudioapi::filesPaneNavigate(file.dir)

# load packages
devtools::load_all("../leaderCluster2")
devtools::load_all("../CNNGP")

# load data
load("GEDI_fitting.Rdata")
load("GEDI_fitting_clustered.Rdata")
load("gedi fitting.Rdata")
load("gedi holdout.Rdata")

data_all = read.csv("./data/me_coords_bio_tcc.csv")
data_missing = data_all[is.na(data_all$bio),]
burnin = 2001

### missing data
class(c.gedi.model) = "NNGP"
clust.missing = predict(c.gedi.model,
                     X.0 = cbind(1,data_missing$tcc),
                     coords.0 = as.matrix(data_missing[,1:2])/1000,
                     sub.sample = list(start = burnin),
                     n.omp.threads = 4,
                     n.report = 1000)
save(clust.missing, file = "cNNGP_pred_missing.Rdata")

class(full.gedi.model) = "NNGP"
full.missing = predict(full.gedi.model,
                     X.0 = cbind(1,data_missing$tcc),
                     coords.0 = as.matrix(data_missing[,1:2])/1000,
                     sub.sample = list(start = burnin),
                     n.omp.threads = 4,
                     n.report = 1000)
save(full.missing, file = "NNNGP_pred_missing.Rdata")

#### held out data
class(c.gedi.model) = "NNGP"
clust.holdout = predict(c.gedi.model,
                        X.0 = cbind(1,holdout.data$tcc),
                        coords.0 = as.matrix(holdout.data[,1:2])/1000,
                        sub.sample = list(start = burnin),
                        n.omp.threads = 4,
                        n.report = 1000)
save(clust.holdout, file = "cNNGP_pred_holdout.Rdata")

class(full.gedi.model) = "NNGP"
full.holdout = predict(full.gedi.model,
                            X.0 = cbind(1,holdout.data$tcc),
                            coords.0 = as.matrix(holdout.data[,1:2])/1000,
                            sub.sample = list(start = burnin),
                            n.omp.threads = 4,
                            n.report = 1000)
save(full.holdout, file = "NNNGP_pred_holdout.Rdata")
