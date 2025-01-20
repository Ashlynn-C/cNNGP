devtools::load_all("../leaderCluster2")
devtools::load_all("../CNNGP")
library(GPvecchia)
library(ggplot2)

for(n in c(2500,10000)){
  
    set.seed(n)
    
    coords = cbind(runif(n,0,1), runif(n,0,1))
    n.keep = 0.8*n
    coords = coords[sample(1:n,n.keep),]
    
    ord = order_maxmin_exact(coords)
    coords_ord = coords[ord,]

    for(m in c(10,20,30)){
      
      indx = mkNNIndxCB(coords_ord, m, 3)
      nnIndx = indx$nnIndx + 1
      nnIndxLU = indx$nnIndxLU
      
      ddNNmat = matrix(NA,ncol=choose((m+1),2),nrow=(n.keep-m))  
      
      for(i in (m+1):n.keep){
        ddm = iDist(coords_ord[c(i,nnIndx[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n.keep+i])]),])
        ddNNmat[(i-m),] = ddm[lower.tri(ddm)]
      }
      
      # get principal components
      pcDDNN = princomp(ddNNmat,cor=T)
      propvar = cumsum(pcDDNN$sdev^2/sum(pcDDNN$sdev^2))
      pc_ddNNmat = pcDDNN$scores[,1:max(1,sum(propvar <= .90))]
      radius_vec_pc = c(0.5,0.75,1:13)
      
      # run clustering for each radius value
      lc_list_pc = lapply(radius_vec_pc, function(x) leaderCluster(pc_ddNNmat, radius = x, max_iter = 1))
      k_vec_pc = sapply(lc_list_pc,'[[', "num_clusters")
      ggplot()+
        geom_point(aes(radius_vec_pc,k_vec_pc))+
        geom_text(aes(radius_vec_pc,k_vec_pc,label=k_vec_pc),hjust=0, vjust=0) + 
        theme_bw() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank()) +
        xlab("Radius") + ylab("Number of clusters") +
        scale_x_continuous(breaks = radius_vec_pc,limits=c(0.4,13.1))

      ggsave(paste0("kappa_sims_plots/","n",n,"_m",m,".png"),width=8,height=5)
  }
}
