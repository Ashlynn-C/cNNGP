# set working directory to source file location
file.dir = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(file.dir))
rstudioapi::filesPaneNavigate(file.dir)

# load packages
devtools::load_all("../leaderCluster2")
devtools::load_all("../CNNGP")
library(GPvecchia)
library(ggplot2)

load("gedi fitting.Rdata")

# set prior on phi
set.seed(98)
n.phi.subsample = 10000
coords = as.matrix(fitting.data[,c("X","Y")])
coords_scaled = coords/1000
dists_scaled = dist(coords_scaled)

min_scaled = min(dists_scaled)
max_scaled = max(dists_scaled)

phi.upper = 3/min_scaled # 2.998
phi.lower = 3/max_scaled # 0.0065

# set model starting places
X = cbind(1,fitting.data$tcc)
Y = fitting.data$bio
beta.ols = solve(t(X)%*%X)%*%t(X)%*%Y
var.est = var(Y-X%*%beta.ols)
tau.sq.starting = 0.5*var.est # give 50% of variance to tau^2
sigma.sq.starting = 0.5*var.est # give 50% of variance to sigma^2
phi.starting = (0.01+3)/2 # mean of prior distribution
starting = list("phi"=phi.starting,"sigma.sq"=sigma.sq.starting,
                "tau.sq"=tau.sq.starting, "beta" = beta.ols)


n.samples = 10000
tuning = list("phi"=0.025)
priors =  list("phi.Unif"=c(0.01,3),"sigma.sq.IG"=c(2, 600), "tau.sq.IG"=c(2,600))
cov.model = "exponential"

# fit clustered model
set.seed(541)
radius = 4
ord = order_maxmin_exact(coords_scaled)
clustering.info = clustering(coords_scaled[ord,], m = 20, radius = radius, omp.threads = 3)
c.gedi.model = clustNNGP(fitting.data$bio~fitting.data$tcc,coords=coords_scaled, ord = ord, radius = radius, 
                                      clust.method = "clustering", PCvar = 0.9, clustering.info = clustering.info, 
                                      starting=starting, method="latent",n.neighbors=20,tuning=tuning,priors=priors,
                                      cov.model=cov.model,n.samples=n.samples, n.omp.threads = 3)

save(c.gedi.model,file="GEDI_fitting_clustered.Rdata")

# fit full model
set.seed(541)
full.gedi.model = clustNNGP(fitting.data$bio~fitting.data$tcc,coords=coords_scaled, ord = ord, clust.method = "exact",
                                     starting=starting, method="latent",n.neighbors=20,tuning=tuning,priors=priors,
                                     cov.model=cov.model,n.samples=n.samples, n.omp.threads = 3)

save(full.gedi.model,file="GEDI_fitting.Rdata")
