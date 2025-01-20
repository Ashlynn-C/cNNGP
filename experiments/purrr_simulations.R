library(purrr)
library(mirai)

# set these values
reps = 30
daemons(4)

# set up simulation scenarios
n.total = c(2500,10000)
n.samples = 50000
burnin = 30000
scenarios.df = expand.grid(n.total = n.total,
                           model = c("NNGP_10", "NNGP_20", "NNGP_30", "cNNGP_10",
                                     "cNNGP_20", "cNNGP_30", "meshedGP", "blockNNGP"),
                          phi = c(-log(0.1)/0.2, -log(0.1)/0.8),
                          rep = 1:reps,
                          stringsAsFactors = F)


radius.df = data.frame(n = rep(n.total, each = 3),
                       m = rep(c(10,30), times = 2),
                       r = c(2,7,2,6))

meshedGP.function = function(data.list, burnin, n.samples){
  
  devtools::load_all("../meshed")
  source("summary_functions.R", local = T)
  
  loc = data.list[["loc"]]
  x.fit = data.list[["x.fit"]]
  y.fit = data.list[["y.fit"]]
  
  # calculate starting points for MCMC models
  keep.idx = (burnin+1):n.samples
  beta.ols = solve(t(x.fit)%*%x.fit)%*%t(x.fit)%*%y.fit
  var.est = var(y.fit-x.fit%*%beta.ols)
  tau.sq.starting = 0.5*var.est # give 50% of variance to tau^2
  sigma.sq.starting = 0.5*var.est # give 50% of variance to sigma^2
  phi.starting = (1+30)/2 # mean of prior distribution
  
  # meshedGP fitting
  set.seed(123)
  meshedGP.time0 = proc.time()
  meshedGP =  spmeshed(y.fit, x.fit, loc,
                       block_size = 30,
                       settings=list(forced_grid = FALSE, ps = FALSE),
                       starting= list(beta = beta.ols, tausq = tau.sq.starting,
                                      theta = matrix(c(phi.starting,sigma.sq.starting),nrow = 2)),
                       prior=list(phi=c(1,30), tausq = c(2,0.1), sigmasq = c(2,1)),
                       n_samples = n.samples,
                       n_burnin = 0,
                       verbose = floor(n.samples/4),
                       n_threads = 1)
  meshedGP.time = proc.time() - meshedGP.time0
  
  # meshedGP prediction
  x.pred = data.list[["x.pred"]]
  y.pred = data.list[["y.pred"]]
  loc.pred = data.list[["loc.pred"]]
  meshedGP.pred.time0 = proc.time()
  meshedGP.pred = predict(meshedGP, x.pred, loc.pred, n_threads = 1, verbose = T)
  meshedGP.pred.time = proc.time() - meshedGP.pred.time0
  
  # order predictions to match holdouts
  meshedGP.pred.order.idx = match(paste(loc.pred[,1], loc.pred[,2]),
                                  paste(meshedGP.pred$coords_out[,1], meshedGP.pred$coords_out[,2]))
  meshedGP.pred.y.samples = meshedGP.pred$preds_out[,1,keep.idx]
  meshedGP.pred.y.samples = meshedGP.pred.y.samples[meshedGP.pred.order.idx,]
  mesehdGP.pred.w.samples = meshedGP.pred$w_out[,keep.idx]
  mesehdGP.pred.w.samples = mesehdGP.pred.w.samples[meshedGP.pred.order.idx,]
  
  # NNGP results
  meshedGP.res.df = results.df(y.pred,
                               meshedGP.pred.y.samples,
                               as.matrix(x.pred%*%meshedGP$beta_mcmc[,1,keep.idx] +
                                           mesehdGP.pred.w.samples),
                               meshedGP$beta_mcmc[1,1,keep.idx] + sapply(meshedGP$w_mcmc[keep.idx], mean),
                               meshedGP$tausq_mcmc[1,keep.idx],
                               meshedGP$lambda_mcmc[1,1,keep.idx]^2,
                               meshedGP$theta_mcmc[1,1,keep.idx] ,
                               meshedGP$beta_mcmc[,1,keep.idx])
  
  meshedGP.res.df["fit.time"] = meshedGP.time["elapsed"]
  meshedGP.res.df["pred.time"] = meshedGP.pred.time["elapsed"]

  return(meshedGP.res.df)
}

blockNNGP.function = function(data.list, keep.idx){
  
  loc = data.list[["loc"]]
  x.fit = data.list[["x.fit"]]
  y.fit = data.list[["y.fit"]]
  w.fit = data.list[["w.fit"]]
  n.partition = 6
  n = nrow(loc)
  source("blockNNGP/blockNNGPfunctionREGULAR.R", local = T)
  source("blockNNGP/blockNNGPrgeneric.R", local = T)
  source("blockNNGP/utils.R", local = T)
  
  library(INLA) 
  library(fields)
  library(lattice)
  library(akima) 
  library(Matrix)
  library(slam)
  library(igraph)
  library(coda)
  library(MBA)
  library(mvtnorm)
  library(ggforce)
  library(Rcpp)
  library(tidyverse)
  library(raster)

  # blockNNGP fitting
  blockNNGP.time0 = proc.time()
  # assign fitting locations to blocks
  blocks.fit = matrix(NA, ncol=1, n)
  hb <- seq(0,1,length.out=n.partition+1)
  vb <- seq(0,1,length.out=n.partition+1)
  k <- 1
  for (j in 2:(n.partition+1)){
    for (i in 2:(n.partition+1)){
      liminfh <- hb[i-1]
      limsuph <- hb[i]
      liminfv <- vb[j-1]
      limsupv <- vb[j]
      indblock <- which(loc[,1] >= liminfh & loc[,1] <= limsuph & loc[,2] >= liminfv & loc[,2] <= limsupv)
      blocks.fit[indblock] <- k
      k <- k + 1
    }
  }
  set.seed(123)
  blockNNGP.model = blockNNGP_reg("regular", loc, y.fit, x.fit, w.fit, blocks.fit, "blockNNGP_plots", n.partition^2, 2)
  blockNNGP.time = proc.time() - blockNNGP.time0
  
  # blockNNGP summary
  x.pred = data.list[["x.pred"]]
  y.pred = data.list[["y.pred"]]
  loc.pred = data.list[["loc.pred"]]
  source("summary_functions.R", local = T)
  blockNNGP.res.df = blockNNGP.summary(blockNNGP.model)
  blockNNGP.res.df["fit.time"] = blockNNGP.time["elapsed"]

  return(blockNNGP.res.df)
}

NNGP.function = function(data.list, m, burnin, n.samples){
  
  library(GPvecchia)
  devtools::load_all("../CNNGP")
  devtools::load_all("../leaderCluster2")
  source("summary_functions.R", local = T)
  
  loc = data.list[["loc"]]
  ord.fit = order_maxmin_exact(loc)
  x.fit = data.list[["x.fit"]]
  y.fit = data.list[["y.fit"]]
  
  # calculate starting points for MCMC models
  keep.idx = (burnin+1):n.samples
  beta.ols = solve(t(x.fit)%*%x.fit)%*%t(x.fit)%*%y.fit
  var.est = var(y.fit-x.fit%*%beta.ols)
  tau.sq.starting = 0.5*var.est # give 50% of variance to tau^2
  sigma.sq.starting = 0.5*var.est # give 50% of variance to sigma^2
  phi.starting = (1+30)/2 # mean of prior distribution
  
  # NNGP fitting
  set.seed(123)
  NNGP.time0 = proc.time()
  NNGP = clustNNGP(y.fit~x.fit-1,coords=loc, ord = ord.fit, clust.method = "exact",
                   starting=list("phi"=phi.starting,"sigma.sq"=sigma.sq.starting,
                                 "tau.sq"=tau.sq.starting, "beta" = beta.ols),
                   method="latent",n.neighbors=m,tuning=list("phi"=0.25),
                   priors=list("phi.Unif"=c(1,30), "sigma.sq.IG"=c(2, 1), "tau.sq.IG"=c(2, 0.1)),
                   cov.model="exponential",n.samples=n.samples, n.omp.threads = 1, n.report = floor(n.samples/5))
  NNGP.time = proc.time() - NNGP.time0
  
  # NNGP prediction
  x.pred = data.list[["x.pred"]]
  y.pred = data.list[["y.pred"]]
  loc.pred = data.list[["loc.pred"]]
  NNGP.pred.time0 = proc.time()
  class(NNGP) = "NNGP" # needed to pass check in prediction code
  NNGP.pred = predict.NNGP(NNGP, x.pred, loc.pred)
  NNGP.pred.time = proc.time() - NNGP.pred.time0
  
  # NNGP results
  NNGP.res.df = results.df(y.pred,
                           NNGP.pred$p.y.0[,keep.idx],
                           x.pred%*%t(as.matrix(NNGP$p.beta.samples)[keep.idx,,drop=FALSE]) +
                             NNGP.pred$p.w.0[,keep.idx,drop=FALSE],
                           NNGP$p.beta.samples[keep.idx,1] + colMeans(NNGP$p.w.samples[,keep.idx]),
                           NNGP$p.theta.samples[keep.idx,"tau.sq"],
                           NNGP$p.theta.samples[keep.idx,"sigma.sq"],
                           NNGP$p.theta.samples[keep.idx,"phi"],
                           t(NNGP$p.beta.samples[keep.idx,]))
  NNGP.res.df["fit.time"] = NNGP.time["elapsed"]
  NNGP.res.df["pred.time"] = NNGP.pred.time["elapsed"]

  return(NNGP.res.df)
}

cNNGP.function = function(data.list, m, radius, burnin, n.samples){
  
  library(GPvecchia)
  devtools::load_all("../CNNGP")
  devtools::load_all("../leaderCluster2")
  source("summary_functions.R", local = T)

  loc = data.list[["loc"]]
  ord.fit = order_maxmin_exact(loc)

  x.fit = data.list[["x.fit"]]
  y.fit = data.list[["y.fit"]]
  
  # calculate starting points for MCMC models
  keep.idx = (burnin+1):n.samples
  beta.ols = solve(t(x.fit)%*%x.fit)%*%t(x.fit)%*%y.fit
  var.est = var(y.fit-x.fit%*%beta.ols)
  tau.sq.starting = 0.5*var.est # give 50% of variance to tau^2
  sigma.sq.starting = 0.5*var.est # give 50% of variance to sigma^2
  phi.starting = (1+30)/2 # mean of prior distribution
  
  # cNNGP fitting
  set.seed(123)
  cNNGP.time0 = proc.time()
  cNNGP = clustNNGP(y.fit~x.fit-1,coords=loc, ord = ord.fit, radius = radius, PC = T, PCvar = 0.9,
                    starting=list("phi"=phi.starting,"sigma.sq"=sigma.sq.starting,
                                  "tau.sq"=tau.sq.starting, "beta" = beta.ols),
                    method="latent",n.neighbors=m,tuning=list("phi"=0.25),
                    priors=list("phi.Unif"=c(1,30), "sigma.sq.IG"=c(2, 1), "tau.sq.IG"=c(2, 0.1)),
                    cov.model="exponential",n.samples=n.samples, n.omp.threads = 1, n.report = floor(n.samples/5))
  cNNGP.time = proc.time() - cNNGP.time0
  
  # cNNGP prediction
  x.pred = data.list[["x.pred"]]
  y.pred = data.list[["y.pred"]]
  loc.pred = data.list[["loc.pred"]]
  cNNGP.pred.time0 = proc.time()
  class(cNNGP) = "NNGP" # needed to pass check in prediction code
  cNNGP.pred = predict.NNGP(cNNGP, x.pred, loc.pred)
  cNNGP.pred.time = proc.time() - cNNGP.pred.time0
  
  # cNNGP results
  cNNGP.res.df = results.df(y.pred,
                            cNNGP.pred$p.y.0[,keep.idx],
                            x.pred%*%t(as.matrix(cNNGP$p.beta.samples)[keep.idx,,drop=FALSE]) +
                              cNNGP.pred$p.w.0[,keep.idx,drop=FALSE],
                            cNNGP$p.beta.samples[keep.idx,1] + colMeans(cNNGP$p.w.samples[,keep.idx]),
                            cNNGP$p.theta.samples[keep.idx,"tau.sq"],
                            cNNGP$p.theta.samples[keep.idx,"sigma.sq"],
                            cNNGP$p.theta.samples[keep.idx,"phi"],
                            t(cNNGP$p.beta.samples[keep.idx,]))
  cNNGP.res.df["fit.time"] = cNNGP.time["elapsed"]
  cNNGP.res.df["pred.time"] = cNNGP.pred.time["elapsed"]

  return(cNNGP.res.df)
}

cat(paste0(Sys.time(), " Started simulations. \n"), file = "simulation_status.txt")
walk(1:nrow(scenarios.df[1:10,]), in_parallel(\(x) {
  
      settings = scenarios.df[x,]
      model = settings[,"model"]
      data.list = readRDS(paste0("datasets/n_",settings["n.total"],
                                 "_phi_", round(settings["phi"],2),
                                 "_rep_", settings["rep"], ".rds"))
      
      # get radius if cNNGP model
      if(endsWith(model, "0") & startsWith(model, "c")) {
        m = substr(model, nchar(model)-1, nchar(model))
        radius = radius.df[which(radius.df$n==settings[,"n.total"] & radius.df$m==m),"r"]
      }
      
      keep.idx = (burnin+1):n.samples
      get.results = function(model) {
          switch(model,
             NNGP_10 = NNGP.function(data.list,10, burnin, n.samples),
             NNGP_20 = NNGP.function(data.list,20, burnin, n.samples),
             NNGP_30 = NNGP.function(data.list,30, burnin, n.samples),
             cNNGP_10 = cNNGP.function(data.list,10,radius, burnin, n.samples),
             cNNGP_20 = cNNGP.function(data.list,20,radius, burnin, n.samples),
             cNNGP_30 = cNNGP.function(data.list,30,radius, burnin, n.samples),
             meshedGP = meshedGP.function(data.list, burnin, n.samples),
             blockNNGP = blockNNGP.function(data.list, keep.idx), 
             stop("Invalid model")
        )
      }
      
      res.df = cbind(data.frame(n = settings[,"n.total"], 
                                phi = round(settings[,"phi"],2),
                                rep = settings[,"rep"],
                                model = model),
                                get.results(model))
      
      try(cat(paste0(Sys.time(), " Completed rep = ", sprintf("%2d", settings[,"rep"]), 
                 ", n = ", sprintf("%5d",settings[,"n.total"]), ", phi = ", sprintf("%5.2f", settings[,"phi"]),
                 ", model = ", sprintf("%9s",model), ".\n"),
          file = "simulation_status.txt", append = T),
          outFile = "try_err.txt")

      save(res.df, file = paste0("model_results/n_",settings[,"n.total"],
                                 "_rep_",settings["rep"], "_phi_", round(settings[,"phi"],2),
                                 "_model_", model, ".Rdata"))
    }, 
    scenarios.df = scenarios.df,
    radius.df = radius.df,
    n.samples = n.samples, 
    burnin = burnin, 
    meshedGP.function = meshedGP.function,
    blockNNGP.function = blockNNGP.function,
    NNGP.function = NNGP.function, 
    cNNGP.function = cNNGP.function
    )
)
daemons(0)
cat(paste0(Sys.time(), " Finished simulations."), file = "simulation_status.txt", append = T)