clustering = function(coords, ord, m, indx, radius = NULL, clust.method = "clustering", PC = T, PCvar = 0.9, max_iter = 1, search.type = NULL, return.object = FALSE, omp.threads = 1){
    

    if(!clust.method %in% c("clustering","exact")){stop("error: Method must be 'clustering' or 'exact'")}
    
    if(!missing(ord)){
        coords = coords[ord,]
    }

    if(missing(indx)){
      ##indexes  
      if(is.null(search.type)){
        indx <- mkNNIndxCB(coords, m, omp.threads)
      }
      if(!is.null(search.type)){
        if(search.type == "brute"){
          indx <- mkNNIndx(coords, m, omp.threads)
        }else{
          indx <- mkNNIndxCB(coords, m, omp.threads)
        }
      }
    }
    
    n = nrow(coords)
    nnIndx = indx$nnIndx + 1
    nnIndxLU = indx$nnIndxLU
    
    ddNNmat = matrix(NA,ncol=choose((m+1),2),nrow=(n-m))  
    
    for(i in (m+1):n){
        ddm = iDist(coords[c(i,nnIndx[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]),])
        ddNNmat[(i-m),]  <- ddm[lower.tri(ddm)]
    }
    
    if(clust.method == "clustering"){
        if(is.null(radius)){stop("error: Radius must be provided if clust.method = clustering")}
        if(PC){
            # get principal components
            pcDDNN = princomp(ddNNmat,cor=T)
            propvar = cumsum(pcDDNN$sdev^2/sum(pcDDNN$sdev^2))
            pc_ddNNmat = pcDDNN$scores[,1:max(1,sum(propvar <= PCvar))]
            lc = leaderCluster(pc_ddNNmat,radius = radius, max_iter = max_iter)
        } else {
            lc = leaderCluster(ddNNmat,radius = radius, max_iter = max_iter)
        }
        
        clustlab = lc$cluster_id
        clust.centers = ddNNmat[lc$cluster_idx,]
        num.clusts = lc$num_clusters
    }
    
    if(clust.method == "exact"){
        clustlab = 1:(n-m)
        clust.centers = ddNNmat
        num.clusts = n-m
    }


    # get distance vectors for point to neighbors
    dist.IU = lapply(1:num.clusts,function(clust){
      dI.U = clust.centers[clust,1:m]
    })
    dist.IU = c(unname(unlist(dist.IU)))

    # get pairwise distance vectors for points in neighbor set
    dist.U = lapply(1:num.clusts,function(clust){
      mmm <- diag(rep(0,(m+1)))
      mmm[lower.tri(mmm)] <- clust.centers[clust,]
      diag(mmm) <- 0
      mmm <- mmm[2:(m+1),2:(m+1)]
      dist.U <- c(mmm + t(mmm))
    })

    dist.U = c(unname(unlist(dist.U)))

    FIndx = c(0:(m-1),clustlab-1+m) # minus 1 for C++ indexing

    BIndx = vector("integer", n)
    BIndx[1] = 0
    for(i in 2:n){
      if(i <= m){
          BIndx[i] = (i-1)*(i-2)/2
       } else {
           BIndx[i] = (m)/2*(m-1)+(clustlab[i-m]-1)*m
       }
    }
    
    if(return.object == TRUE){
        return(list(dist.IU = dist.IU, dist.U = dist.U, clustIndx = clustlab-1, BIndx = BIndx, FIndx = FIndx, k = num.clusts, lc = lc))
    } else{
        return(list(dist.IU = dist.IU, dist.U = dist.U, clustIndx = clustlab-1, BIndx = BIndx, FIndx = FIndx, k = num.clusts))
    }
}
