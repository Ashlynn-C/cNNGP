library(purrr)
library(mirai)

# set these two values
daemons(4)
reps = 30

# data settings
n.vec = c(2500, 10000)
phi.vec = c(-log(0.1)/0.2, -log(0.1)/0.4, -log(0.1)/0.8)
sigma.sq = 1
B = as.matrix(c(1,5))
tau.sq = 0.1

rmvn = function(n, mu=0, V = matrix(1)){
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
        stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

data.seeds = data.frame(n = rep(n.vec, each = 3*reps),
                        phi = rep(phi.vec, each = reps, times = 2),
                        rep = rep(1:reps, times = 6),
                        seed = 1:(reps*6))


for(n.total in n.vec){
  for(phi in phi.vec) {
    walk(1:reps, in_parallel(\(rep) {

      # generate data 
      set.seed(data.seeds[(data.seeds$n == n.total) & (data.seeds$phi == phi) & (data.seeds$rep == rep), "seed"])
      coords = cbind(runif(n.total,0,1), runif(n.total,0,1))
      x = cbind(1, rnorm(n.total))
      D = as.matrix(dist(coords))
      R = exp(-phi*D)
      w = rmvn(1, rep(0,n.total), sigma.sq*R)
      y = rnorm(n.total, x%*%B + w, sqrt(tau.sq))
      
      # partition data into fitting/testing with polygons
      block.data = data.frame(x = coords[,1], y = coords[,2], value = y)
      block.data.sf = sf::st_as_sf(block.data, coords = c("x", "y"), crs = "EPSG:3857")
      folds = blockCV::cv_spatial(x = block.data.sf,
                         k = 5,
                         rows_cols = c(20, 20),
                         selection = "random", #
                         iteration = 1, 
                         biomod2 = F,
                         plot = F,
                         report = F) 
      
      fit.idx = folds$folds_list[[1]][[1]]
      
      loc = coords[fit.idx,]
      x.fit = x[fit.idx,]
      y.fit = y[fit.idx]
      w.fit = w[fit.idx]

      loc.pred = coords[-fit.idx,]
      x.pred = x[-fit.idx,]
      y.pred = y[-fit.idx]
      w.pred = w[-fit.idx]
      
      data.list = list(loc = loc, x.fit = x.fit, y.fit = y.fit, w.fit = w.fit, 
                       loc.pred = loc.pred, x.pred = x.pred, y.pred = y.pred, w.pred = w.pred)
      
      saveRDS(data.list, 
           file = paste0("datasets/n_",n.total,"_phi_", round(phi,2),"_rep_", rep, ".rds"))
    },
    rmvn = rmvn,
    n.total = n.total, 
    phi = phi,
    sigma.sq = sigma.sq, 
    B = B, 
    tau.sq = tau.sq,
    data.seeds = data.seeds
    ))
  }
}

daemons(0)
