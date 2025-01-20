# set working directory to source file location
file.dir = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(file.dir))
# load packages
devtools::load_all("leaderCluster2")
devtools::load_all("CNNGP")

##Make some data
set.seed(3)
n <- 300
coords <- cbind(runif(n,0,1), runif(n,0,1))

x <- cbind(1, rnorm(n))
B <- as.matrix(c(1,5))

sigma.sq <- 5
tau.sq <- 1
phi <- 3/0.5

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, x%*%B + w, sqrt(tau.sq))

##Fit a clustered and sequential NNGP model
n.samples <- 500

starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)

tuning <- list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5)

priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))

cov.model <- "exponential"

m = 20
# clustNNGP function with clust.method = "exact" is the same as the spNNGP function
set.seed(3)
model.clust = clustNNGP(y~x-1,coords=coords,starting=starting, method="latent",n.neighbors=m, clust.method = "exact",
                         tuning=tuning,priors=priors,cov.model=cov.model,n.samples=n.samples, n.omp.threads = 3)

set.seed(3)
model.full = spNNGP(y~x-1,coords=coords,starting=starting,method="latent",n.neighbors=m, 
                        tuning=tuning,priors=priors,cov.model=cov.model,n.samples=n.samples, n.omp.threads = 3)


plot(model.clust$p.theta.samples)
model.clust$run.time
plot(model.full$p.theta.samples)
model.full$run.time