# set working directory to source file location
file.dir = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(file.dir))
rstudioapi::filesPaneNavigate(file.dir)

# load packages
devtools::load_all("../leaderCluster2")
devtools::load_all("../CNNGP")
library(GPvecchia)
library(ggplot2)
library(dplyr)

load("data/gedi.training.Rdata")
coords = as.matrix(gedi.training[,c("x","y")])
coords_scaled = coords/1000
ord = order_maxmin_exact(coords_scaled)
coords_ord = coords_scaled[ord,]

# make sure to match search type to that used in CNNGP fitting
m = 20
indx = mkNNIndxCB(coords_ord, m, 2)

n = nrow(coords_ord)
nnIndx = indx$nnIndx + 1
nnIndxLU = indx$nnIndxLU

ddNNmat = matrix(NA,ncol=choose((m+1),2),nrow=(n-m))  

for(i in (m+1):n){
    ddm = iDist(coords_ord[c(i,nnIndx[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]),])
    ddNNmat[(i-m),]  <- ddm[lower.tri(ddm)]
}

# get principal components
pcDDNN = princomp(ddNNmat,cor=T)
propvar = cumsum(pcDDNN$sdev^2/sum(pcDDNN$sdev^2))
pc_ddNNmat = pcDDNN$scores[,1:max(1,sum(propvar <= .9))] # 35 columns

# set number of points to subsample for clustering 
set.seed(123)
num.subsample = 10000
sub_idx = sample(1:nrow(pc_ddNNmat), num.subsample, replace = F)

dists_sub = pc_ddNNmat[sub_idx,]
radius_vec = 1:10

# run clustering for each radius value
lc_list = lapply(radius_vec, function(x) leaderCluster(dists_sub, radius = x, max_iter = 1))

# plot result
k_vec = sapply(lc_list,'[[', "num_clusters")
ggplot()+geom_point(aes(radius_vec,k_vec))+    
    geom_text(aes(radius_vec,k_vec,label=k_vec),hjust=-0.15, vjust=0) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +
    xlab("Radius") + ylab("Number of clusters") + 
    scale_x_continuous(breaks = seq(1,10,1),limits=c(0.95,10.5)) +
    scale_y_continuous(breaks = c(0,2000,4000,6000),limits=c(100,7100))
ggsave("gedi_kappa_PC.png",width=5,height=3)
